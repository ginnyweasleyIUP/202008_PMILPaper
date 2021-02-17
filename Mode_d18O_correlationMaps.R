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
source("~/R_function/base_fun.R")

## functions
##############

load_data <- function(exp) {
  # Isotopes
  print("Isotopes")
  file <- paste0("/modeldata/hadcm3/isotopes/monthly_fixed/", exp, ".nc")
  ncfile <- nc_open(file)
  tryCatch(lon <<- ncvar_get(ncfile, "longitude"), error = function(e) {lon <<- ncvar_get(ncfile, "longitude_1")})
  tryCatch(lat <<- ncvar_get(ncfile, "latitude"), error = function(e) {lat <<- ncvar_get(ncfile, "latitude_1")})
  time <- ncvar_get(ncfile, "t")
  data <- ncvar_get(ncfile, "dO18")
  nc_close(ncfile)
  data_list <- list(lon = lon, lat= lat, time = time, dO18 = data)
  
  # Precipitation
  print("Precipitation")
  file <- paste0("/modeldata/hadcm3/precipitation/monthly_fixed/", exp, ".nc")
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
  file <- paste0("/modeldata/hadcm3/surface_temperature/monthly_fixed/", exp, ".nc")
  ncfile <- nc_open(file)
  lon <- ncvar_get(ncfile, "longitude")
  lat <- ncvar_get(ncfile, "latitude")
  time <- ncvar_get(ncfile, "t")
  data <- ncvar_get(ncfile, "temp_1")
  nc_close(ncfile)
  data_list$temp <- data
  
  # sea level pressure
  print("pressure")
  file <- paste0("/modeldata/hadcm3/sea_level_pressure/monthly_fixed/", exp, ".nc")
  ncfile <- nc_open(file)
  lon <- ncvar_get(ncfile, "longitude")
  lat <- ncvar_get(ncfile, "latitude")
  time <- ncvar_get(ncfile, "t")
  data <- ncvar_get(ncfile, "p")
  nc_close(ncfile)
  sel <- which(data_list$time %in% time)
  data_list$psl <- array(NA, dim=c(dim(data)[1], dim(data)[2], length(data_list$time)))
  data_list$psl[,,sel] <- data
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


plot_cor <- function(data, breaks) {
  
  p <- ggplot(data = data, aes(x, y)) + coord_quickmap() +
    geom_raster(aes(fill=pmin(pmax(cor, min(breaks)), max(breaks))), interpolate=T) +
    #geom_polygon(data = Map(0, 360, -90, 90, 8), aes(x = long, y = lat, group = group), 
    #             colour = "black", fill = "transparent", size = .5) +
    #borders("world", wrap = c(0, 360), ylim = c(-90, 90)) +
    geom_polygon(data = map_data('world', wrap=c(0,360), ylim=c(-90,90)), aes(x=long, y = lat, group = group), 
                 fill="transparent", colour = "black") +
    geom_contour(aes(x, y, z=cor), breaks=breaks, size=.4, linetype = "22", colour = "black") +
    #theme + theme1 +
    theme(panel.background = element_rect(fill = "transparent", colour= "black"),
          plot.margin = margin(15,15,15,15),
          plot.title   = element_text(size = 24, face = "bold",
                                      margin = margin(0,0,20,0), hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(face = "bold", size=24),
          legend.text = element_text(size = 20),

          axis.text    = element_text(size = 14),
          axis.title   = element_text(size = 18),
          axis.title.x = element_text(margin = margin(20,0,0,0)),
          axis.title.y = element_text(margin = margin(0,20,0,0))) +
    scale_fill_gradientn(colours =  rev(brewer.pal(name="RdBu", n=11)), limits=range(breaks), breaks = breaks,
                         name = "Correlation", na.value="grey",
                         guide = guide_colourbar(title.position = "top", title.hjust=.5, title.vjust=0,  order = 1,
                                                 barwidth = 20, barheight = 1.5)) +
    labs(x=NULL, y=NULL) +
    scale_x_continuous(expand = c(0,0), limits = c(0, 360), breaks = seq(0, 360, 30)) +
    scale_y_continuous(expand = c(0,0), limits = c(-90, 90), breaks = seq(-90, 90, 30))
  
  return(p)
}


plot_enso <- function(data_list, exp, t_start, t_end, variable) {
  print("nino Index")
  sellon_Nino3.4 <- which(data_list$lon >= 360-170 & data_list$lon <= 360-120)
  sellat_Nino3.4 <- which(data_list$lat >= -5 & data_list$lat <= 5)
  Nino3.4 <- apply(data_list$temp[sellon_Nino3.4, sellat_Nino3.4, ], 3, mean)
  Nino3.4_anom <- anomm(Nino3.4)
  Nino3.4_anom_DJF <- rollapply(Nino3.4_anom[(t_start+11):t_end], 3, mean, by=12)
  #Nino3.4_anom_SONDJF <- rollapply(Nino3.4_anom[1:12000+8], 6, mean, by=12)
  
  p <- list()
  
  
  if("temperature" %in% variable) {
    print("temperature")
    data_plot <- data.frame(cor = c(apply(apply(data_list$temp[,,(t_start+11):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = Nino3.4_anom_DJF, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_temp <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="surface temperature anomaly VS Nino3.4 (DJF) ")
    p<- c(p, list(p_temp))
  }
  
  if("precipitation" %in% variable) {
    print("precipitation")
    data_plot <- data.frame(cor = c(apply(apply(data_list$pr[,,(t_start+11):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = Nino3.4_anom_DJF, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_pr <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="precipitation anomaly VS Nino3.4 (DJF) ")
    p<- c(p, list(p_pr))
  }
  
  if("isotope" %in% variable) {
    print("isotope")
    data_plot <- data.frame(cor = c(apply(apply(data_list$dO18[,,(t_start+11):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = Nino3.4_anom_DJF, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_dO18 <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="d18O anomaly VS Nino3.4 (DJF) ")
    p<- c(p, list(p_dO18))
  }


  p_all <- ggarrange(plotlist=p, nrow=length(variable), ncol=1, common.legend = T, legend = "bottom")
  #ggsave(paste0("~/08_heidelberg_2020-2022/Corals/Figures/ENSO_cor_", exp, ".pdf"), p, device="pdf", width=10, height=20)
  return(p_all)
}



plot_ISM <- function(data_list, exp, t_start, t_end, variable) {
  print("South Asian Monsoon Index")
  sellon_ISM <- which(data_list$lon >= 60 & data_list$lon <= 100)
  sellat_ISM <- which(data_list$lat >= 10 & data_list$lat <= 30)
  ISM <- apply(data_list$pr[sellon_ISM, sellat_ISM, ], 3, mean)
  ISM_anom <- anomm(ISM)
  ISM_anom_JJAS <- rollapply(ISM_anom[(t_start+5):t_end], 4, mean, by=12)
  
  p <- list()
  
  
  if("temperature" %in% variable) {
    print("temperature")
    data_plot <- data.frame(cor = c(apply(apply(data_list$temp[,,(t_start+5):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = ISM_anom_JJAS, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_temp <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="surface temperature anomaly VS ISM (JJAS) ")
    p<- c(p, list(p_temp))
  }
  
  if("precipitation" %in% variable) {
    print("precipitation")
    data_plot <- data.frame(cor = c(apply(apply(data_list$pr[,,(t_start+5):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = ISM_anom_JJAS, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_pr <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="precipitation anomaly VS ISM (JJAS) ")
    p<- c(p, list(p_pr))
  }
  
  if("isotope" %in% variable) {
    print("isotope")
    data_plot <- data.frame(cor = c(apply(apply(data_list$dO18[,,(t_start+5):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = ISM_anom_JJAS, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_dO18 <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="d18O anomaly VS ISM (JJAS) ")
    p<- c(p, list(p_dO18))
  }
  
  
  p_all <- ggarrange(plotlist=p, nrow=length(variable), ncol=1, common.legend = T, legend = "bottom")
  #ggsave(paste0("~/08_heidelberg_2020-2022/Corals/Figures/ISM_cor_", exp, ".pdf"), p, device="pdf", width=10, height=20)
  return(p_all)
}


plot_NAO <- function(data_list, exp, t_start, t_end, variable) {
  print("NAO")
  sellon_NAO <- which(data_list$lon >= 360-90 | data_list$lon <= 40)
  sellat_NAO <- which(data_list$lat >= 20 & data_list$lat <= 80)
  selnoNA_DJF <- unlist(sapply(1:trunc((t_end-t_start-1)/12), function(i) {
    sel <- 12:14+(i-1)*12
    if (all(!is.na(data_list$psl[1,1,sel]))) return(sel)
  }))
  
  EOF_DJF <- prcomp(apply(t(matrix(data_list$psl[sellon_NAO,sellat_NAO,selnoNA_DJF], length(sellon_NAO)*length(sellat_NAO), length(selnoNA_DJF))),
                          2, anomm), scale=T, center=T); gc()
#  plot the first PC
#  test = data.frame(cor = EOF_DJF$rotation[,1],
#    x   = rep(data_list$lon[sellon_NAO], length(data_list$lat[sellat_NAO])),
#    y   = rep(data_list$lat[sellat_NAO], each=length(data_list$lon[sellon_NAO]))            
#  )
  
  
#  plot_cor(test, seq(-1,1,.2)*0.05)
  
  NAO_DJF <- EOF_DJF$x[,1]
  NAO_DJF <- NAO_DJF*sign(cor(NAO_DJF, data_list$psl[which(data_list$lon == 345), which(data_list$lat == 65), selnoNA_DJF]))
  
  NAO_mDJF <- rollapply(NAO_DJF, 3, mean, by=3)
  
  p <- list()
  
  
  if("temperature" %in% variable) {
    print("temperature")
    data_plot <- data.frame(cor = c(apply(apply(data_list$temp[,,(t_start+11):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = NAO_mDJF, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_temp <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="surface temperature anomaly VS NAO (DJF) ")
    p<- c(p, list(p_temp))
  }
  
  if("precipitation" %in% variable) {
    print("precipitation")
    data_plot <- data.frame(cor = c(apply(apply(data_list$pr[,,(t_start+11):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = NAO_mDJF, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_pr <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="precipitation anomaly VS NAO (DJF) ")
    p<- c(p, list(p_pr))
  }
  
  if("pressure" %in% variable) {
    print("presure")
    data_plot <- data.frame(cor = c(apply(apply(data_list$psl[,,(t_start+11):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = NAO_mDJF, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_psl <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="psl anomaly VS NAO (DJF) ")
    p<- c(p, list(p_psl))
  }
  
  if("isotope" %in% variable) {
    print("isotope")
    data_plot <- data.frame(cor = c(apply(apply(data_list$dO18[,,(t_start+11):t_end], 1:2, function(x)
      rollapply(anomm(x), 3, mean, by=12)), 2:3, cor, y = NAO_mDJF, use = "na")),
      x   = rep(data_list$lon, length(data_list$lat)),
      y   = rep(data_list$lat, each=length(data_list$lon)))
    p_dO18 <- plot_cor(data_plot, seq(-1,1,.2)) + labs(title="d18O anomaly VS NAO (DJF) ")
    p<- c(p, list(p_dO18))
  }
  
  
  p_all <- ggarrange(plotlist=p, nrow=2, ncol=1, common.legend = T, legend = "bottom")
  #ggsave(paste0("~/08_heidelberg_2020-2022/Corals/Figures/NAO_cor_", exp, ".pdf"), p, device="pdf", width=10, height=20)
  return(p_all)
}


################################################################################
################################################################################

## Run functions for each simulation
######################################



data_list <- load_data("xnapa")
p_ENSO <- plot_enso(data_list, "xnapa", 1, 12000, c("temperature", "isotope"))
p_ISM  <- plot_ISM(data_list, "xnapa", 1, 12000, c("precipitation", "isotope"))
p_NAO  <- plot_NAO(data_list, "xnapa", 1, 12000, c("pressure", "isotope"))




data_list <- load_data("xmzke")
plot_enso(data_list, "xmzke", 12, 12011)
data_list <- load_data("xmzkj")
plot_enso(data_list, "xmzkj", 12, 12011)
data_list <- load_data("xnage")
plot_enso(data_list, "xnage", 12, 12011)
data_list <- load_data("xnapb")
plot_enso(data_list, "xnapb", 12, 12011)
data_list <- load_data("xnhca")
plot_enso(data_list, "xnhca", 12, 12011)
data_list <- load_data("xnuaa")
plot_enso(data_list, "xnuaa", 12, 12011)
data_list <- load_data("xnuav")
plot_enso(data_list, "xnuav", 12, 599)
data_list <- load_data("xmzka")
plot_enso(data_list, "xmzka", 12, 12011)
data_list <- load_data("xmzkf")
plot_enso(data_list, "xmzkf", 12, 263)
data_list <- load_data("xmzkk")
plot_enso(data_list, "xmzkk", 12, 12011)
data_list <- load_data("xnagf")
plot_enso(data_list, "xnagf", 12, 12011)
data_list <- load_data("xnapc")
plot_enso(data_list, "xnapc", 12, 12011)
data_list <- load_data("xnhcb")
plot_enso(data_list, "xnhcb", 12, 11219)
data_list <- load_data("xnuab")
plot_enso(data_list, "xnuab", 12, 539)
data_list <- load_data("xnuca")
plot_enso(data_list, "xnuca", 12, 599)
data_list <- load_data("xmzkb")
plot_enso(data_list, "xmzkb", 12, 12011)
data_list <- load_data("xmzkg")
plot_enso(data_list, "xmzkg", 12, 12011)
data_list <- load_data("xnagb")
plot_enso(data_list, "xnagb", 12, 12011)
data_list <- load_data("xnagg")
plot_enso(data_list, "xnagg", 12, 12011)
data_list <- load_data("xnapo")
plot_enso(data_list, "xnapo", 12, 12011)
data_list <- load_data("xnhdb")
plot_enso(data_list, "xnhdb", 12, 12011)
data_list <- load_data("xnuac")
plot_enso(data_list, "xnuac", 12, 599)
data_list <- load_data("xnucb")
plot_enso(data_list, "xnucb", 12, 671)
data_list <- load_data("xmzkc")
plot_enso(data_list, "xmzkc", 12, 12011)
data_list <- load_data("xmzkh")
plot_enso(data_list, "xmzkh", 12, 12011)
data_list <- load_data("xnagc")
plot_enso(data_list, "xnagc", 12, 2759)
data_list <- load_data("xnagh")
plot_enso(data_list, "xnagh", 12, 12011)
data_list <- load_data("xnapp")
plot_enso(data_list, "xnapp", 12, 12011)
data_list <- load_data("xnhdc")
plot_enso(data_list, "xnhdc", 12, 10055)
data_list <- load_data("xnuad")
plot_enso(data_list, "xnuad", 12, 12011)
data_list <- load_data("xnucc")
plot_enso(data_list, "xnucc", 12, 599)
data_list <- load_data("xmzkd")
plot_enso(data_list, "xmzkd", 12, 12011)
data_list <- load_data("xmzki")
plot_enso(data_list, "xmzki", 12, 599)
data_list <- load_data("xnagd")
plot_enso(data_list, "xnagd", 12, 12011)
data_list <- load_data("xnapq")
plot_enso(data_list, "xnapq", 12, 12011)
data_list <- load_data("xnhdd")
plot_enso(data_list, "xnhdd", 12, 2603)
data_list <- load_data("xnuae")
plot_enso(data_list, "xnuae", 12, 599)
data_list <- load_data("xnuaa")
plot_enso(data_list, "xnuaa", 12, 599)
