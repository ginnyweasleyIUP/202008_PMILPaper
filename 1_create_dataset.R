#################################################
## 1 Create Data Set ############################

library(plyr)
library(dplyr)
library(stacy.hadcm.tools)
library(PaleoSpec)
library(nest)
library(tidyverse)

LOCAL = TRUE
d.lon = 96
d.lat = 73

#################################################

DATA_past1000 <- list()
for(run in c("a","b","c")){
  DATA_past1000[[paste0("SIM_yearly_",run)]] <- list(
    TEMP = list(),
    PREC = list(), 
    ISOT = list()
  )
}

DATA_past1000$CAVES <- list()
if(LOCAL){
  DATA_past1000$CAVES$site_info <- read.csv("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/02_SISAL/SISAL_v2/site_countries.csv")
}else{
  DATA_past1000$CAVES$site_info <- read.csv("/stacywork/ginnyweasley/02_SISAL/SISAL_v2/site_countries.csv")  
}

DATA_past1000$CAVES$entity_info <- list()
DATA_past1000$CAVES$record_data <- list()
DATA_past1000$CAVES$sim_data_yearly <- list()
DATA_past1000$CAVES$sim_data_downsampled <- list()


DATA_past1000_SIM_RAW <- list()
for(run in c("a","b","c")){
  DATA_past1000_SIM_RAW[[paste0(run)]] <- list(
    TEMP = list(),
    PREC = list(), 
    ISOT = list()
  )
}

DATA_past1000$SIM_mean <- list()

#################################################

source("Functions/clear_data_matrix.R")

for(run in c("a", "b", "c")){
  print(run)
  print("TEMP")
  if(LOCAL){
    ncf <- (ncdf4::nc_open(paste0("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnap", run, "_surface_temperature.nc")))
  }else{
    ncf <- (ncdf4::nc_open(paste0("/stacywork/hadcm3/surface_temperature/monthly_fixed/xnap", run, ".nc"))) 
  }
  DATA_past1000_SIM_RAW[[run]]$TEMP <- clear_data_matrix(ncdf4::ncvar_get(ncf),1)
  DATA_past1000_SIM_RAW[[run]]$TEMP_t <- ncf$dim$t$vals
  ncdf4::nc_close(ncf)
  
  print("PREC")
  if(LOCAL){
    ncf<-ncdf4::nc_open(paste0("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnap",run,"_precipitation.nc")) 
  }else{
    ncf<-ncdf4::nc_open(paste0("/stacywork/hadcm3/precipitation/monthly_fixed/xnap",run,".nc")) 
  }
  DATA_past1000_SIM_RAW[[run]]$PREC <- clear_data_matrix(ncdf4::ncvar_get(ncf),2)
  DATA_past1000_SIM_RAW[[run]]$PREC_t <- ncf$dim$t$vals
  DATA_past1000$SIM_mean$lon <- ncf$dim$longitude$vals
  DATA_past1000$SIM_mean$lat <- ncf$dim$latitude$vals
  ncdf4::nc_close(ncf)
  print("ISOT")
  if(LOCAL){
    ncf<-ncdf4::nc_open(paste0("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnap",run,"_isotopes.nc"))
  }else{
    ncf<-ncdf4::nc_open(paste0("/stacywork/hadcm3/isotopes/monthly_fixed/xnap",run,".nc"))
  }
  DATA_past1000_SIM_RAW[[run]]$ISOT <- clear_data_matrix(ncdf4::ncvar_get(ncf, 'dO18'),3)
  DATA_past1000_SIM_RAW[[run]]$ISOT_t <- ncf$dim$t$vals
  ncdf4::nc_close(ncf)
  
  if(LOCAL){
    ncf<-ncdf4::nc_open(paste0("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/05_HadCM3/xnap/xnap",run,"_sea_level_pressure.nc"))
  }else{
    ncf<-ncdf4::nc_open(paste0("/stacywork/hadcm3/sea_level_pressure/monthly_fixed/xnap",run,".nc"))
  }
  DATA_past1000_SIM_RAW[[run]]$SLPR <- ncdf4::ncvar_get(ncf)
  DATA_past1000_SIM_RAW[[run]]$SLPR_t <- ncf$dim$t$vals
  ncdf4::nc_close(ncf)
  
}


if(LOCAL){
  ncf<-ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/06_Daten/05_HadCM3/hadcm_oro_001cyr.nc")
}else{
  ncf<-ncdf4::nc_open("/stacywork/hadcm3/hadcm_oro_001cyr.nc")
}

DATA_past1000$SIM_mean$OROG <- ncdf4::ncvar_get(ncf)
ncdf4::nc_close(ncf)

remove(ncf)

remove(run)
remove(clear_data_matrix)

#################################################
## Align timeseries #############################
#################################################

year_start = 850
year_stop = 1850

for(run in c("a","b","c")){
  for(var in c("TEMP", "PREC", "ISOT", "SLPR")){
    #we want to start in March
    pos_start = length(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]][DaysSinceToAD(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]])<year_start]) -1+4
    pos_stop = length(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]]) - length(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]][DaysSinceToAD(DATA_past1000_SIM_RAW[[run]][[paste0(var, "_t")]])>year_stop]) -1+4
    
    DATA_past1000_SIM_RAW[[run]][[var]] <- DATA_past1000_SIM_RAW[[run]][[var]][,,pos_start:pos_stop]
    DATA_past1000_SIM_RAW[[run]][[paste0(var,"_t")]] <- DATA_past1000_SIM_RAW[[run]][[paste0(var,"_t")]][pos_start:pos_stop]
  }
  
}

rm(run,var, pos_start, pos_stop)

DATA_past1000$time <- c(year_start, year_stop)
remove(year_start, year_stop)


#################################################
## Raw Caves

source("Functions/SISAL_extracting.R")

data <- load_sisal_data(year_start = DATA_past1000$time[1], year_stop = DATA_past1000$time[2])
DATA_past1000$CAVES$entity_info <- data[[1]] %>% group_by(entity_id) %>% arrange(entity_id)

#Schmeißt alle Höhlen raus, die nicht gebraucht werden
DATA_past1000$CAVES$site_info <- DATA_past1000$CAVES$site_info %>% filter(site_id %in% DATA_past1000$CAVES$entity_info$site_id)
for (ii in DATA_past1000$CAVES$entity_info$entity_id){
  name = paste0("ENTITY", ii)
  #if(ii%%10 == 0){
  print(name)
  #}
  site <- DATA_past1000$CAVES$entity_info %>% filter(entity_id == ii) %>% distinct(site_id)
  DATA_past1000$CAVES$record_data[[name]] <- data[[2]] %>% filter(entity_id == ii) %>% distinct(entity_id, mineralogy, interp_age, d18O_measurement) %>%
    mutate(site_id = (site$site_id))
  #if(ii == 144 | ii == 226){next}
  #DATA_past1000$CAVES$record_data[[name]]$chron <- as.tibble(data[[3]] %>% filter(entity_id == ii))
}

remove(data, site, ii, name, load_sisal_data)

#################################################
## Extract Cave Data 

source("Functions/extract_gridboxes.R")

for(run in c("a", "b", "c")){
  print(run)
  for (ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
    lon_cave = DATA_past1000$CAVES$entity_info$longitude[ii]
    
    if(lon_cave<0){lon_cave = 360+lon_cave}
    
    lat_cave = DATA_past1000$CAVES$entity_info$latitude[ii]
    entity_id = DATA_past1000$CAVES$entity_info$entity_id[ii]
    
    ratios <- extract_gridboxes(lon_cave, lat_cave)
    
    name <- paste0("ENTITY",entity_id)
    DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_",run)]] <- rowSums(cbind(ratios$E1*DATA_past1000_SIM_RAW[[run]]$TEMP[ratios$E1_lon_pos, ratios$E1_lat_pos,],
                                                                                     ratios$E2*DATA_past1000_SIM_RAW[[run]]$TEMP[ratios$E2_lon_pos, ratios$E2_lat_pos,],
                                                                                     ratios$E3*DATA_past1000_SIM_RAW[[run]]$TEMP[ratios$E3_lon_pos, ratios$E3_lat_pos,],
                                                                                     ratios$E4*DATA_past1000_SIM_RAW[[run]]$TEMP[ratios$E4_lon_pos, ratios$E4_lat_pos,]), na.rm = T)
    
    DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_",run)]] <- rowSums(cbind(ratios$E1*DATA_past1000_SIM_RAW[[run]]$PREC[ratios$E1_lon_pos, ratios$E1_lat_pos,],
                                                                                     ratios$E2*DATA_past1000_SIM_RAW[[run]]$PREC[ratios$E2_lon_pos, ratios$E2_lat_pos,],
                                                                                     ratios$E3*DATA_past1000_SIM_RAW[[run]]$PREC[ratios$E3_lon_pos, ratios$E3_lat_pos,],
                                                                                     ratios$E4*DATA_past1000_SIM_RAW[[run]]$PREC[ratios$E4_lon_pos, ratios$E4_lat_pos,]), na.rm = T)
    
    DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_",run)]] <- rowSums(cbind(ratios$E1*DATA_past1000_SIM_RAW[[run]]$ISOT[ratios$E1_lon_pos, ratios$E1_lat_pos,],
                                                                                     ratios$E2*DATA_past1000_SIM_RAW[[run]]$ISOT[ratios$E2_lon_pos, ratios$E2_lat_pos,],
                                                                                     ratios$E3*DATA_past1000_SIM_RAW[[run]]$ISOT[ratios$E3_lon_pos, ratios$E3_lat_pos,],
                                                                                     ratios$E4*DATA_past1000_SIM_RAW[[run]]$ISOT[ratios$E4_lon_pos, ratios$E4_lat_pos,]), na.rm = T)
  }
}


elevation_cave_sim <- array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id), 3))
colnames(elevation_cave_sim) <- c("entity_id", "elevation_sim", "sim-cave")

for (ii in 1:(length(DATA_past1000$CAVES$entity_info$entity_id))){
  lon_cave = DATA_past1000$CAVES$entity_info$longitude[ii]
  if(lon_cave<0){lon_cave = 360+lon_cave}
  lat_cave = DATA_past1000$CAVES$entity_info$latitude[ii]
  
  ratios <- extract_gridboxes(lon_cave, lat_cave)
  
  elevation_cave_sim[ii,1] <- DATA_past1000$CAVES$entity_info$entity_id[ii]
  elevation_cave_sim[ii,2] <- ratios$E1*DATA_past1000$SIM_mean$OROG[ratios$E1_lon_pos, ratios$E1_lat_pos] +
    ratios$E2*DATA_past1000$SIM_mean$OROG[ratios$E2_lon_pos, ratios$E2_lat_pos] +
    ratios$E3*DATA_past1000$SIM_mean$OROG[ratios$E3_lon_pos, ratios$E3_lat_pos] +
    ratios$E4*DATA_past1000$SIM_mean$OROG[ratios$E4_lon_pos, ratios$E4_lat_pos]
  
  elevation_cave_sim[ii,3] <- elevation_cave_sim[ii,2] - DATA_past1000$CAVES$entity_info$elevation[ii]
  
  
}

DATA_past1000$CAVES$entity_info <- right_join(DATA_past1000$CAVES$entity_info, as.data.frame(elevation_cave_sim), by = 'entity_id')

remove(ratios, ii, lat_cave, lon_cave, name, entity_id, elevation_cave_sim, run, extract_gridboxes)

################################################
## Annual data

for(run in c("a","b","c")){
  DATA_past1000[[paste0("SIM_yearly_",run)]]$TEMP <- array(dim = c(d.lon,d.lat,diff(DATA_past1000$time)))
  DATA_past1000[[paste0("SIM_yearly_",run)]]$PREC <- array(dim = c(d.lon,d.lat,diff(DATA_past1000$time)))
  DATA_past1000[[paste0("SIM_yearly_",run)]]$ISOT <- array(dim = c(d.lon,d.lat,diff(DATA_past1000$time)))
  DATA_past1000[[paste0("SIM_yearly_",run)]]$ITPC <- array(dim = c(d.lon,d.lat,diff(DATA_past1000$time))) # <- prec weighted mean
  
  for (lon in (1:d.lon)){
    for (lat in 1:d.lat){
      print(paste(lon,lat))
      for(year in 1:diff(DATA_past1000$time)){
        pos_start = 12*(year-1)+1
        pos_stop  = 12*(year-1)+12
        DATA_past1000[[paste0("SIM_yearly_",run)]]$TEMP[lon,lat,year] = mean(DATA_past1000_SIM_RAW[[run]]$TEMP[lon,lat, pos_start:pos_stop], na.rm = T)
        DATA_past1000[[paste0("SIM_yearly_",run)]]$PREC[lon,lat,year] = mean(DATA_past1000_SIM_RAW[[run]]$PREC[lon,lat, pos_start:pos_stop], na.rm = T)
        DATA_past1000[[paste0("SIM_yearly_",run)]]$ISOT[lon,lat,year] = mean(DATA_past1000_SIM_RAW[[run]]$ISOT[lon,lat, pos_start:pos_stop], na.rm = T)
        DATA_past1000[[paste0("SIM_yearly_",run)]]$ITPC[lon,lat,year] = sum(DATA_past1000_SIM_RAW[[run]]$PREC[lon,lat, pos_start:pos_stop]*DATA_past1000_SIM_RAW[[run]]$ISOT[lon,lat, pos_start:pos_stop],
                                                                            na.rm = T)/sum(DATA_past1000_SIM_RAW[[run]]$PREC[lon,lat, pos_start:pos_stop], na.rm = T)
      }
    }
  }
}


remove(lon,lat,year, pos_start, pos_stop, run)

for(run in c("a", "b", "c")){
  print(run)
  for (ii in 1:(length(DATA_past1000$CAVES$entity_info$entity_id))){
    print(ii)
    entity_id = DATA_past1000$CAVES$entity_info$entity_id[ii]
    name <- paste0("ENTITY", entity_id)
    DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("TEMP_", run)]] <- numeric(diff(DATA_past1000$time))
    DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("PREC_", run)]] <- numeric(diff(DATA_past1000$time))
    DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ISOT_", run)]] <- numeric(diff(DATA_past1000$time))
    DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ITPC_", run)]] <- numeric(diff(DATA_past1000$time))
    for(year in 1:diff(DATA_past1000$time)){
      pos_start = 12*(year-1)+1
      pos_stop  = 12*(year-1)+12
      DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("TEMP_", run)]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop], na.rm = T)
      DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("PREC_", run)]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop], na.rm = T)
      DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ISOT_", run)]][year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop], na.rm = T)
      DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ITPC_", run)]][year] <- sum(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*
                                                                                         DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop], na.rm = T)/
        sum(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop], na.rm = T)
    }
  }
}


remove(entity_id, year, pos_start, pos_stop, name, ii, run)


#################################################
## Down-sample Cave Data

source("Functions/SubsampleTimeseriesBlock_highresNA.R")

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  print(ii)
  name = paste0("ENTITY", DATA_past1000$CAVES$entity_info$entity_id[ii])
  for(run in c("a", "b", "c")){
    assign(paste0("data_temp_", run), SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("TEMP_",run)]]),
                                                                            start = 1950-DATA_past1000$time[2],
                                                                            end = 1950-DATA_past1000$time[1]),
                                                                         DATA_past1000$CAVES$record_data[[name]]$interp_age))
    
    assign(paste0("data_prec_", run), SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("PREC_",run)]]),
                                                                            start = 1950-DATA_past1000$time[2],
                                                                            end = 1950-DATA_past1000$time[1]),
                                                                         DATA_past1000$CAVES$record_data[[name]]$interp_age))
    
    assign(paste0("data_isot_", run), SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ISOT_",run)]]),
                                                                            start = 1950-DATA_past1000$time[2],
                                                                            end = 1950-DATA_past1000$time[1]),
                                                                         DATA_past1000$CAVES$record_data[[name]]$interp_age))
    
    assign(paste0("data_itpc_", run), SubsampleTimeseriesBlock_highresNA(ts(data = rev(DATA_past1000$CAVES$sim_data_yearly[[name]][[paste0("ITPC_",run)]]),
                                                                            start = 1950-DATA_past1000$time[2],
                                                                            end = 1950-DATA_past1000$time[1]),
                                                                         DATA_past1000$CAVES$record_data[[name]]$interp_age))
  }
  data <- matrix(c(DATA_past1000$CAVES$record_data[[name]]$interp_age, 
                   data_temp_a, data_prec_a, data_isot_a, data_itpc_a,
                   data_temp_b, data_prec_b, data_isot_b, data_itpc_b,
                   data_temp_c, data_prec_c, data_isot_c, data_itpc_c), ncol= 13)
  colnames(data) = c("interp_age", 
                     "TEMP_a", "PREC_a","ISOT_a", "ITPC_a", 
                     "TEMP_b", "PREC_b","ISOT_b", "ITPC_b", 
                     "TEMP_c", "PREC_c","ISOT_c", "ITPC_c")
  
  DATA_past1000$CAVES$sim_data_downsampled[[name]] <- as.tibble(data)
}

remove(name, ii, data,
       data_isot_a, data_itpc_a, data_prec_a, data_temp_a,
       data_isot_b, data_itpc_b, data_prec_b, data_temp_b,
       data_isot_c, data_itpc_c, data_prec_c, data_temp_c, run)
remove(SubsampleTimeseriesBlock_highresNA)

#################################################
## Seasonal Cave Data 

# winter_mask = c( 1, 1, 1,NA,NA,NA,NA,NA,NA,NA,NA,NA)
# spring_mask = c(NA,NA,NA, 1, 1, 1,NA,NA,NA,NA,NA,NA)
# summer_mask = c(NA,NA,NA,NA,NA,NA, 1, 1, 1,NA,NA,NA)
# autumn_mask = c(NA,NA,NA,NA,NA,NA,NA,NA,NA, 1, 1, 1)


spring_mask = c( 1, 1, 1,NA,NA,NA,NA,NA,NA,NA,NA,NA)
summer_mask = c(NA,NA,NA, 1, 1, 1,NA,NA,NA,NA,NA,NA)
autumn_mask = c(NA,NA,NA,NA,NA,NA, 1, 1, 1,NA,NA,NA)
winter_mask = c(NA,NA,NA,NA,NA,NA,NA,NA,NA, 1, 1, 1)

#DATA_past1000$CAVES$sim_data_seasonal <- vector(mode = "list")

for(run in c("a", "b", "c")){
  for (ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
    entity_id = DATA_past1000$CAVES$entity_info$entity_id[ii]
    name = paste0("ENTITY", entity_id)
    
    DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]] = list(
      SUMMER = list(temp_mean = numeric(diff(DATA_past1000$time)), prec_mean = numeric(diff(DATA_past1000$time)), isot_mean = numeric(diff(DATA_past1000$time))),
      AUTUMN = list(temp_mean = numeric(diff(DATA_past1000$time)), prec_mean = numeric(diff(DATA_past1000$time)), isot_mean = numeric(diff(DATA_past1000$time))),
      WINTER = list(temp_mean = numeric(diff(DATA_past1000$time)), prec_mean = numeric(diff(DATA_past1000$time)), isot_mean = numeric(diff(DATA_past1000$time))),
      SPRING = list(temp_mean = numeric(diff(DATA_past1000$time)), prec_mean = numeric(diff(DATA_past1000$time)), isot_mean = numeric(diff(DATA_past1000$time)))
    )
    
    
    for(year in 1:diff(DATA_past1000$time)){
      pos_start = 12*(year-1)+1
      pos_stop  = 12*(year-1)+12
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$SUMMER$temp_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop]*summer_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$SUMMER$prec_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop]*summer_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$SUMMER$isot_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*summer_mask, na.rm = T)
      
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$AUTUMN$temp_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop]*autumn_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$AUTUMN$prec_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop]*autumn_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$AUTUMN$isot_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*autumn_mask, na.rm = T)
      
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$WINTER$temp_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop]*winter_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$WINTER$prec_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop]*winter_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$WINTER$isot_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*winter_mask, na.rm = T)
      
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$SPRING$temp_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("TEMP_", run)]][pos_start:pos_stop]*spring_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$SPRING$prec_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("PREC_", run)]][pos_start:pos_stop]*spring_mask, na.rm = T)
      DATA_past1000$CAVES$sim_data_seasonal[[run]][[name]]$SPRING$isot_mean[year] <- mean(DATA_past1000$CAVES$sim_data_raw[[name]][[paste0("ISOT_", run)]][pos_start:pos_stop]*spring_mask, na.rm = T)
      
    }
  }
}

remove(ii, name, pos_start, pos_stop, year, entity_id, run)

#################################################
## ARAGONITE, CALCITE, SMOW 

## Aragonite and Calcite d18O have to be converted to drip water equivalents to be comparable
## Also in ORder to be comparable to SMOW which is the standard of the simulation we need to convert the dripwater d18O from VPDB to SMOW
print("Drip-Water Conversion")

mineralogy = character(length(DATA_past1000$CAVES$entity_info$entity_id))

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]

  print(entity)
  data_rec <- DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]
  data_sim <- DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]
  dw_eq_a <- numeric(length(data_sim$interp_age))
  dw_eq_b <- numeric(length(data_sim$interp_age))
  dw_eq_c <- numeric(length(data_sim$interp_age))
  for(run in c("a", "b", "c")){
    dw_eq <- numeric(length(data_sim$interp_age))
    for(jj in 1:length(data_sim$interp_age)){
      if(data_rec$mineralogy[jj] == "calcite"){
        dw_eq[jj] = 1.03092 * (data_rec$d18O_measurement[jj] - ((16.1*1000)/(data_sim[[paste0("TEMP_", run)]][jj]+273.15)-24.6)) + 30.92
      }else if(data_rec$mineralogy[jj] == "aragonite"){
        dw_eq[jj] = 1.03092 * (data_rec$d18O_measurement[jj] - ((18.34*1000)/(data_sim[[paste0("TEMP_", run)]][jj]+273.15)-31.954)) + 30.92
      }else{
        dw_eq[jj] = NA
      }
    }
    assign(paste0("dw_eq_",run), dw_eq)
  }
  
  min_test = data_rec %>% group_by(mineralogy) %>% count()
  if(dim(min_test)[1] == 1){
    mineralogy[ii] = min_test$mineralogy[1]
  }else{
    mineralogy[ii] = "mixed"
  }
  
  mineralogy[ii] 
  
  DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]] <- data.frame(
    site_id = data_rec$site_id,
    entity_id = data_rec$entity_id,
    interp_age = data_rec$interp_age,
    d18O_measurement = data_rec$d18O_measurement,
    d18O_dw_eq_a = dw_eq_a,
    d18O_dw_eq_b = dw_eq_b,
    d18O_dw_eq_c = dw_eq_c
  )
}

DATA_past1000$CAVES$entity_info$mineralogy <- mineralogy


remove(dw_eq_a, dw_eq_b, dw_eq_c, ii, data_rec, data_sim, dw_eq, jj, entity)

#################################################
## UMWANDLUNGEN

DATA_past1000$CAVES$site_info <- DATA_past1000$CAVES$site_info %>% mutate(elevation = as.numeric(as.character(elevation)))

#################################################
## Timeseries

load("Data/Timeseries.RData")

source("Functions/aw_mean.R")
value <- list("a" = numeric(1000), "b" = numeric(1000), "c" = numeric(1000))
for(ii in 1:length(value$a)){
  value$a[ii] <- simpleawmean(DATA_past1000$SIM_yearly_a$TEMP[,,ii])
  value$b[ii] <- simpleawmean(DATA_past1000$SIM_yearly_b$TEMP[,,ii])
  value$c[ii] <- simpleawmean(DATA_past1000$SIM_yearly_c$TEMP[,,ii])
}
#mean(value$a)*500000/(4419*mean(value$a)-500000) Temp Änderung pro Delta Änderung

Timeseries$HadCM3_TAM <- NULL
Timeseries$HadCM3_GMST_a <- zoo(x = value$a-mean(value$a[1:100])-0.2, order.by = seq(from = -1*(1140), to = -1*(100), by = 1))
Timeseries$HadCM3_GMST_b <- zoo(x = value$b-mean(value$b[1:100])-0.2, order.by = seq(from = -1*(1140), to = -1*(100), by = 1))
Timeseries$HadCM3_GMST_c <- zoo(x = value$c-mean(value$c[1:100])-0.2, order.by = seq(from = -1*(1140), to = -1*(100), by = 1))

##Bunker site_id = 117, entity_id = 240,242
# HadCM3 Bunker cave --> site 117
data_240 <- DATA_past1000$CAVES$record_data$ENTITY240
data_242 <- DATA_past1000$CAVES$record_data$ENTITY242

Timeseries$SISAL_Bunker_240 <- zoo(x = data_240$d18O_dw_eq_a, order.by = -1*data_240$interp_age)
Timeseries$SISAL_Bunker_242 <- zoo(x = data_242$d18O_dw_eq_a, order.by = -1*data_242$interp_age)
data_240 <- DATA_past1000$CAVES$sim_data_yearly$ENTITY240
Timeseries$HadCM3_Bunker_a <- zoo(x = data_240$ISOT_a, order.by = seq(from = -1*(1100)+1, to = -1*(100), by = 1)) 
Timeseries$HadCM3_Bunker_b <- zoo(x = data_240$ISOT_b, order.by = seq(from = -1*(1100)+1, to = -1*(100), by = 1))
Timeseries$HadCM3_Bunker_c <- zoo(x = data_240$ISOT_c, order.by = seq(from = -1*(1100)+1, to = -1*(100), by = 1))

rm(value, data_240, data_242, ii)

save(Timeseries, file = "Data/Timeseries.RData")


#################################################
## DATA EXPORT

# Yearly Data

for(run in c("a", "b", "c")){
  DATA_EXPORT_YEARLY <- list()
  for(entity in DATA_past1000$CAVES$entity_info$entity_id){
    data_new = array(dim = c(length(DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$TEMP_a),7))
    colnames(data_new) = c("site_id", "entity_id", "year_BP", "TEMP", "PREC", "ISOT", "ITPC")
    data_new[,1] = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
    data_new[,2] = entity
    data_new[,3] = seq(1950-DATA_past1000$time[1], 1950-DATA_past1000$time[2]+1, by = -1)
    data_new[,4] = DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$TEMP_a
    data_new[,5] = DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$PREC_a
    data_new[,6] = DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$ISOT_a
    data_new[,7] = DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$ITPC_a
    
    if(entity == 14){data = data_new}
    else{data = rbind(data, data_new)}
    
  }
  DATA_EXPORT_YEARLY <- data
  
  write.csv(DATA_EXPORT_YEARLY, file = paste0("Data/SISAL_HadCM3_xnap",run,"_PMIL_yearly.csv"), row.names = F)
  
}

# Seasonal Data
for(run in c("a", "b", "c")){
  print(run)
  DATA_EXPORT_SEASONAL <- list()
  for(season in c("WINTER", "SPRING", "SUMMER", "AUTUMN")){
    print(season)
    for(entity in DATA_past1000$CAVES$entity_info$entity_id){
      print(entity)
      data_new = array(dim = c(length(DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$TEMP_a),6))
      colnames(data_new) = c("site_id", "entity_id", "year_BP", "TEMP", "PREC", "ISOT")
      data_new[,1] = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
      data_new[,2] = entity
      data_new[,3] = seq(1950-DATA_past1000$time[1], 1950-DATA_past1000$time[2]+1, by = -1)
      data_new[,4] = DATA_past1000$CAVES$sim_data_seasonal[[run]][[paste0("ENTITY", entity)]][[season]]$temp_mean
      data_new[,5] = DATA_past1000$CAVES$sim_data_seasonal[[run]][[paste0("ENTITY", entity)]][[season]]$prec_mean
      data_new[,6] = DATA_past1000$CAVES$sim_data_seasonal[[run]][[paste0("ENTITY", entity)]][[season]]$isot_mean

      if(entity == 14){data = data_new}
      else{data = rbind(data, data_new)}
      
    }
    DATA_EXPORT_SEASONAL <- data
    print("here")
    write.csv(DATA_EXPORT_SEASONAL, file = paste0("Data/Seasonal/SISAL_HadCM3_xnap",run,"_PMIL_",season,".csv"), row.names = F) 
  }
}


# Down-Sampled Data

DATA_EXPORT_DS <- list()
for(entity in DATA_past1000$CAVES$entity_info$entity_id){
  data_new = array(dim = c(length(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age),19))
  colnames(data_new) = c("site_id", "entity_id", "interp_age", "d18O_measurement", 
                         "d18O_dw_eq_a","d18O_dw_eq_b","d18O_dw_eq_c", 
                         "TEMP_a", "PREC_a", "ISOT_a", "ITPC_a", "TEMP_b", "PREC_b", "ISOT_b", "ITPC_b", "TEMP_c", "PREC_c", "ISOT_c", "ITPC_c")
  data_new[,1] = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  data_new[,2] = entity
  data_new[,3] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$interp_age
  data_new[,4] = DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_measurement
  
  data_new[,5] = DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_a
  data_new[,6] = DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_b
  data_new[,7] = DATA_past1000$CAVES$record_data[[paste0("ENTITY", entity)]]$d18O_dw_eq_c
  
  data_new[,8] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP_a
  data_new[,9] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC_a
  data_new[,10] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT_a
  data_new[,11] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_a

  data_new[,12] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP_b
  data_new[,13] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC_b
  data_new[,14] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT_b
  data_new[,15] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_b

  data_new[,16] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP_c
  data_new[,17] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC_c
  data_new[,18] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT_c
  data_new[,19] = DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_c
  
  if(entity == 14){data = data_new}
  else{data = rbind(data, data_new)}
  
}
DATA_EXPORT_DS <- data

write.csv(DATA_EXPORT_DS, file = paste0("Data/SISAL_HadCM3_ds.csv"), row.names = F) 

# Entity Info

DATA_EXPORT_INFO <- list()
data = array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id),19))
colnames(data) = c("site_id", "entity_id", "latitude", "longitude", "elevation", "geology", "cover_thickness", "distance_entrance", "elevation_sim", "mineralogy",
                   "temp_ds", "prec_ds", "isot_ds", "itps_ds", "temp_full", "prec_full", "isot_full", "itps_full", "winter_prec")
data[,1] = DATA_past1000$CAVES$entity_info$site_id
data[,2] = DATA_past1000$CAVES$entity_info$entity_id
data[,3] = DATA_past1000$CAVES$entity_info$latitude
data[,4] = DATA_past1000$CAVES$entity_info$longitude
data[,5] = DATA_past1000$CAVES$entity_info$elevation
data[,6] = DATA_past1000$CAVES$entity_info$geology
data[,7] = DATA_past1000$CAVES$entity_info$cover_thickness
data[,8] = DATA_past1000$CAVES$entity_info$distance_entrance
data[,9] = round(DATA_past1000$CAVES$entity_info$elevation_sim, digits = 4)
data[,10] = DATA_past1000$CAVES$entity_info$mineralogy
for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  data[ii,11] = round(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP_a, na.rm = T), digits = 4)
  data[ii,12] = round(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC_a*8.6148e4, na.rm = T), digits = 4)
  data[ii,13] = round(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT_a, na.rm = T), digits = 4)
  data[ii,14] = round(mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_a, na.rm = T), digits = 4)
  
  data[ii,15] = round(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$TEMP_a, na.rm = T), digits = 4)
  data[ii,16] = round(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$PREC_a*8.6148e4, na.rm = T), digits = 4)
  data[ii,17] = round(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$ISOT_a, na.rm = T), digits = 4)
  data[ii,18] = round(mean(DATA_past1000$CAVES$sim_data_yearly[[paste0("ENTITY", entity)]]$ITPC_a, na.rm = T), digits = 4)
  
  data[ii,19] = round(mean(DATA_past1000$CAVES$sim_data_seasonal$a[[paste0("ENTITY", entity)]]$WINTER$prec_mean*8.6148e4, na.rm = T), digits = 4)
  
  # data[ii,10] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP_b, na.rm = T)
  # data[ii,10] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC_b, na.rm = T)
  # data[ii,10] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT_b, na.rm = T)
  # data[ii,10] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_b, na.rm = T)
  # data[ii,10] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$TEMP_c, na.rm = T)
  # data[ii,10] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$PREC_c, na.rm = T)
  # data[ii,10] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ISOT_c, na.rm = T)
  # data[ii,10] = mean(DATA_past1000$CAVES$sim_data_downsampled[[paste0("ENTITY", entity)]]$ITPC_c, na.rm = T)
  

  
}
DATA_EXPORT_INFO <- data

write.csv(DATA_EXPORT_INFO, file = paste0("Data/SISAL_HadCM3_entity_info.csv"), row.names = F) 


# yearly HadCM3 --> for Mean Value and for Correlation

#################################################
## already export fields for plotting...

DATA_EXPORT_FIELD <- list("a" <- list(),
                          "b" <- list(),
                          "c" <- list())
for(run in c("a", "b", "c")){
  DATA_EXPORT_FIELD[[run]]$TEMP_MEAN <- apply(DATA_past1000[[paste0("SIM_yearly_",run)]]$TEMP, c(1,2), mean, na.rm = T)
  DATA_EXPORT_FIELD[[run]]$PREC_MEAN <- apply(DATA_past1000[[paste0("SIM_yearly_",run)]]$PREC, c(1,2), mean, na.rm = T)
  DATA_EXPORT_FIELD[[run]]$ISOT_MEAN <- apply(DATA_past1000[[paste0("SIM_yearly_",run)]]$ISOT, c(1,2), mean, na.rm = T)
  DATA_EXPORT_FIELD[[run]]$ITPC_MEAN <- apply(DATA_past1000[[paste0("SIM_yearly_",run)]]$ITPC, c(1,2), mean, na.rm = T)
}

save(DATA_EXPORT_FIELD, file = "Data/LM_HadCM3_annualmean.RData")

DATA_EXPORT_CORR <- list("a" <- list(),
                         "b" <- list(),
                         "c" <- list())

for(run in c("a", "b", "c")){
  CORR_FIELD <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)),
                     TEMP_P = array(dim = c(96,73)), PREC_P = array(dim = c(96,73)))
  for(var in c("TEMP", "PREC")){
    for(lon in 1:96){
      print(lon)
      for(lat in 1:73){
        if(sum(is.na(DATA_past1000[[paste0("SIM_yearly_", run)]]$ISOT[lon,lat,]))<100){
          CORR = cor.test(DATA_past1000[[paste0("SIM_yearly_", run)]][[var]][lon,lat,], DATA_past1000[[paste0("SIM_yearly_", run)]]$ITPC[lon,lat,], na.rm = TRUE)
          CORR_FIELD[[var]][lon,lat] = CORR$estimate[[1]]
          CORR_FIELD[[paste0(var, "_P")]][lon,lat] = CORR$p.value
        }else{
          CORR_FIELD[[var]][lon,lat] = NA
          CORR_FIELD[[paste0(var, "_P")]][lon,lat] = NA
        }
      }
    }
  }
  DATA_EXPORT_CORR[[run]] <- CORR_FIELD
}

save(DATA_EXPORT_CORR, file = "Data/LM_HadCM3_correlation.RData")

# VEGETATION

DATA_VEGETATION <- list()
nc <- ncdf4::nc_open("/home/ginnyweasley/Dokumente/01_Promotion/07_R_Code/202005_Paper/Data/HadCM3/xnapa_pi_0_decade_0.nc")
lon <- ncdf4::ncvar_get(nc,"longitude")
lat <- ncdf4::ncvar_get(nc,"latitude")
pseudo <- ncdf4::ncvar_get(nc,"pseudo")
surface <- ncdf4::ncvar_get(nc,"surface")
time <- ncdf4::ncvar_get(nc,"t")
pftfrac <- ncdf4::ncvar_get(nc,"field1391")
pft_max <- apply(pftfrac[,,1:9,6],c(1,2),function(x) {if (!is.na(x[1])) {return(which.max(x))} else {NA}})
pft_max[is.na(pft_max)] = 7 #turn into water
DATA_VEGETATION <- pft_max

save(DATA_VEGETATION, file = "Data/LM_HadCM3_JuneVegetation.RData")

rm(d.lat, d.lon, entity, ii, run, season)
rm(data, DATA_EXPORT_DS, DATA_EXPORT_INFO, DATA_EXPORT_SEASONAL, DATA_EXPORT_YEARLY, data_new)
rm(a,b,c,CORR, CORR_FIELD, DATA_EXPORT_CORR, DATA_EXPORT_FIELD, DATA_VEGETATION, min_test, nc, pft_max)
rm(lat, LOCAL, lon, mineralogy, pftfrac, pseudo, spring_mask, summer_mask, surface, time)
rm(autumn_mask, var, winter_mask)
