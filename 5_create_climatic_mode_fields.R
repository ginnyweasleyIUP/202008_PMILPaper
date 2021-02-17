################################################################################
################################################################################
####################### Modes of variability in hadcm3 #########################
################################################################################
################################################################################

## libraries
##############
require(ggplot2)
require(ncdf4)
require(zoo)
require(RColorBrewer)
require(ggpubr)
library(stacy.hadcm.tools)
source("Functions/base_fun.R")

## functions
##############

load_data <- function(exp) {
  # Isotopes
  print("Isotopes")
  file <- paste0("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/",exp,"_isotopes.nc")
  ncfile <- nc_open(file)
  tryCatch(lon <<- ncvar_get(ncfile, "longitude"), error = function(e) {lon <<- ncvar_get(ncfile, "longitude_1")})
  tryCatch(lat <<- ncvar_get(ncfile, "latitude"), error = function(e) {lat <<- ncvar_get(ncfile, "latitude_1")})
  time <- ncvar_get(ncfile, "t")
  data <- ncvar_get(ncfile, "dO18")
  nc_close(ncfile)
  data_list <- list(lon = lon, lat= lat, time = time, dO18 = data)
  
  # Precipitation
  print("Precipitation")
  file <- paste0("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/",exp,"_precipitation.nc")
  ncfile <- nc_open(file)
  lon <- ncvar_get(ncfile, "longitude")
  lat <- ncvar_get(ncfile, "latitude")
  time <- ncvar_get(ncfile, "t")
  data <- ncvar_get(ncfile, "precip")
  nc_close(ncfile)
  #if (all.equal(data_list$lat,lat) & all.equal(data_list$lon,lon) & all.equal(data_list$time,time))
  data_list$pr <- data
  
  # surface temperature (including sst)
  print("temperature")
  file <- paste0("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/",exp,"_surface_temperature.nc")
  ncfile <- nc_open(file)
  lon <- ncvar_get(ncfile, "longitude")
  lat <- ncvar_get(ncfile, "latitude")
  time <- ncvar_get(ncfile, "t")
  data <- ncvar_get(ncfile, "temp_1")
  nc_close(ncfile)
  data_list$temp <- data
  
  # sea level pressure
  print("pressure")
  file <- paste0("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnapa_sea_level_pressure.nc")
  ncfile <- nc_open(file)
  lon <- ncvar_get(ncfile, "longitude")
  lat <- ncvar_get(ncfile, "latitude")
  time <- ncvar_get(ncfile, "t")
  data <- ncvar_get(ncfile, "p")
  nc_close(ncfile)
  sel <- which(data_list$time %in% time)
  data_list$psl <- array(NA, dim=c(dim(data)[1], dim(data)[2], length(data_list$time)))
  data_list$psl[,,sel] <- data
  
  pos = which(data_list$time == (ADToDaysSince(850)+15))
  data_list$time <- data_list$time[pos:(pos+12000-1)]
  data_list$temp <- data_list$temp[,,pos:(pos+12000-1)]
  data_list$pr   <- data_list$pr[,,pos:(pos+12000-1)]
  data_list$psl  <- data_list$psl[,,pos:(pos+12000-1)]
  data_list$dO18 <- data_list$dO18[,,pos:(pos+12000-1)]
  
  rm(pos,sel,data,time,lon,lat,ncfile,file)
  
  return(data_list)
}

anomm <- function(x) { #remove seasonal mean from monthly data
  y <- x
  for (m in 1:12) {
    selm <- seq(m, length(y), 12)
    y[selm] <- y[selm] - mean(y[selm], na.rm=T)
  }
  return(y)
}

ENSO_field <- function(data_list, t_start, t_end){
  #d18O anomaly VS Nino3.4 (DJF)
  sellon_Nino3.4 <- which(data_list$lon >= 360-170 & data_list$lon <= 360-120)
  sellat_Nino3.4 <- which(data_list$lat >= -5 & data_list$lat <= 5)
  Nino3.4 <- apply(data_list$temp[sellon_Nino3.4, sellat_Nino3.4, ], 3, mean)
  Nino3.4_anom <- anomm(Nino3.4)
  Nino3.4_anom_DJF <- rollapply(Nino3.4_anom[(t_start+11):t_end], 3, mean, by=12)
  
  CORR <- array(dim = c(96,73))+NA
  P <- array(dim = c(96,73))+NA
  for(lon in 1:96){
    for(lat in 2:72){
      c <- cor.test(rollapply(data_list$dO18[lon,lat,(t_start+11):(t_end)],3,mean, by = 12), Nino3.4_anom_DJF, na.rm = T)
      CORR[lon,lat] = c$estimate[[1]]
      P[lon,lat] = c$p.value
    }
  }
  
  results = list(corr = CORR, p = P)
  
  return(results)
}

ISM_field <- function(data_list, t_start, t_end){
  #d18O anomaly VS ISM (JJAS)
  sellon_ISM <- which(data_list$lon >= 60 & data_list$lon <= 100)
  sellat_ISM <- which(data_list$lat >= 10 & data_list$lat <= 30)
  ISM <- apply(data_list$pr[sellon_ISM, sellat_ISM, ], 3, mean)
  ISM_anom <- anomm(ISM)
  ISM_anom_JJAS <- rollapply(ISM_anom[(t_start+5):t_end], 4, mean, by=12)
  
  CORR <- array(dim = c(96,73))+NA
  P <- array(dim = c(96,73))+NA
  for(lon in 1:96){
    for(lat in 2:72){
      c <- cor.test(rollapply(data_list$dO18[lon,lat,(t_start+5):(t_end)],4,mean, by = 12), ISM_anom_JJAS, na.rm = T)
      CORR[lon,lat] = c$estimate[[1]]
      P[lon,lat] = c$p.value
    }
  }
  
  results = list(corr = CORR, p = P)
  
  return(results)
  
}

NAO_field <- function(data_list, t_start, t_end){
  #d18O anomaly VS NAO (DJF)
  sellon_NAO <- which(data_list$lon >= 360-90 | data_list$lon <= 40)
  sellat_NAO <- which(data_list$lat >= 20 & data_list$lat <= 80)
  selnoNA_DJF <- unlist(sapply(1:trunc((t_end-t_start-1)/12), function(i) {
    sel <- 12:14+(i-1)*12
    if (all(!is.na(data_list$psl[1,1,sel]))) return(sel)
  }))
  
  EOF_DJF <- prcomp(apply(t(matrix(data_list$psl[sellon_NAO,sellat_NAO,selnoNA_DJF], length(sellon_NAO)*length(sellat_NAO), length(selnoNA_DJF))),
                          2, anomm), scale=T, center=T); gc()
  
  NAO_DJF <- EOF_DJF$x[,1]
  NAO_DJF <- NAO_DJF*sign(cor(NAO_DJF, data_list$psl[which(data_list$lon == 345), which(data_list$lat == 65), selnoNA_DJF]))
  
  NAO_mDJF <- rollapply(NAO_DJF, 3, mean, by=3)
  
  CORR <- array(dim = c(96,73))+NA
  P <- array(dim = c(96,73))+NA
  for(lon in 1:96){
    for(lat in 2:72){
      c <- cor.test(rollapply(data_list$dO18[lon,lat,(t_start+11):(t_end)],3,mean, by = 12)[1:994], NAO_mDJF, na.rm = T)
      CORR[lon,lat] = c$estimate[[1]]
      P[lon,lat] = c$p.value
    }
  }
  results = list(corr = CORR, p = P)
  
  return(results)
  
}

####MAIN#############################################################
experiment = "xnapa"

ClimModes <- list(ENSO = list(), NAO = list(), ISM = list())
data_list <- load_data(experiment)
ClimModes$ENSO <- ENSO_field(data_list, 1, 12000)
ClimModes$NAO  <- NAO_field(data_list, 1, 12000)
ClimModes$ISM  <- ISM_field(data_list, 1, 12000)

save(ClimModes, file = paste0("ClimModes_",experiment,".RData"))

