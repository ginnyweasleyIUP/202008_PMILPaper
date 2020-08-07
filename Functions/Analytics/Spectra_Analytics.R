#################################################
## ANALYSIS TO SPECTRA 

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

SPECTRA <- list(
  RECORDS = list(),
  SIM_ds = list("a" = list(),"b" = list(),"c" = list()),
  SIM_full = list("a" = list(),"b" = list(),"c" = list()),
  MEAN_SPEC = list()
)

# spectrum needs equally spaced data...


# Proxies
entities_spec <- list()
entities_lats <- list()

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])){
  print(ii)
  entity = DATA_past1000$CAVES$entity_info$entity_id[mask_spec][ii]
  data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
  if(length(data_rec$interp_age) > 8 & mask_spec[ii] & entity != 351){
    #85 -> eID351 der Quatsch macht allgemein
    entities_spec = c(entities_spec, entity)
    entities_lats = c(entities_lats, DATA_past1000$CAVES$entity_info$latitude[ii])
    start_ts = ceiling(head(data_rec$interp_age, n = 1))
    length = length(data_rec$interp_age)
    stop_ts = floor(tail(data_rec$interp_age, n = 1))
    if(length > (stop_ts-start_ts)){length = (stop_ts-start_ts)}
    stop_ts = floor((stop_ts-start_ts)/length)*length+start_ts
    record <- PaleoSpec::MakeEquidistant(data_rec$interp_age,data_rec$d18O_measurement,time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data <- ts(data = record, start = start_ts, end = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    name = paste0("ENTITY", entity)
    
    SPECTRA$RECORDS[[name]] = spectrum(data, plot = F)
    
  }
  
}

SPECTRA$entities_spec_rec <- as.numeric(entities_spec)

SPECTRA$MEAN_SPEC[["Record"]] <- MeanSpectrum(SPECTRA$RECORDS)
SPECTRA$MEAN_SPEC_WEIGH[["Record"]] <- MeanSpectrum(SPECTRA$RECORDS, weights = cos(as.numeric(entities_lats)*pi/180)/sum(cos(as.numeric(entities_lats)*pi/180)))

# Down-sampled

for(run in c("a","b","c")){
  SPECTRA$SIM_ds[[run]] <- list(
    ISOT = list(), 
    ITPC = list()
  )
  
  for(ii in 1:length(SPECTRA$entities_spec_rec)){
    print(ii)
    entity = SPECTRA$entities_spec_rec[ii]
    name = paste0("ENTITY", entity)
    data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
    length = length(data_rec$interp_age)
    start_ts = ceiling(head(data_rec$interp_age, n = 1))
    stop_ts = floor(tail(data_rec$interp_age, n = 1))
    if(length(data_rec$interp_age) > (stop_ts-start_ts)){length = (stop_ts-start_ts)}
    
    stop_ts = floor((stop_ts-start_ts)/length)*length+start_ts
    
    #ISOT
    record <- PaleoSpec::MakeEquidistant(data_rec$interp_age,data_rec[[paste0("ISOT_", run)]], time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    SPECTRA$SIM_ds[[run]]$ISOT[[name]] = spectrum(data, plot = F)
    
    #ITPC
    record <- PaleoSpec::MakeEquidistant(data_rec$interp_age,data_rec[[paste0("ITPC_", run)]],time.target = seq(from = start_ts, to = stop_ts, by = floor((stop_ts-start_ts)/length)))
    record = na.omit(record)
    data  = ts(data = record, start = start_ts, end   = stop_ts, deltat = floor((stop_ts-start_ts)/length))
    
    SPECTRA$SIM_ds[[run]]$ITPC[[name]] = spectrum(data, plot = F)
    
  }
  
  SPECTRA$MEAN_SPEC$ISOT$ds[[paste0("SIM_ds_ISOT_", run)]] <- MeanSpectrum(SPECTRA$SIM_ds[[run]]$ISOT)
  SPECTRA$MEAN_SPEC$ITPC$ds[[paste0("SIM_ds_ITPC_", run)]] <- MeanSpectrum(SPECTRA$SIM_ds[[run]]$ITPC)
  
  SPECTRA$MEAN_SPEC_WEIGH$ISOT$ds[[paste0("SIM_ds_ISOT_", run)]] <- MeanSpectrum(SPECTRA$SIM_ds[[run]]$ISOT, 
                                                                                 weights = cos(as.numeric(entities_lats)*pi/180)/sum(cos(as.numeric(entities_lats)*pi/180)))
  SPECTRA$MEAN_SPEC_WEIGH$ITPC$ds[[paste0("SIM_ds_ITPC_", run)]] <- MeanSpectrum(SPECTRA$SIM_ds[[run]]$ITPC, 
                                                                                 weights = cos(as.numeric(entities_lats)*pi/180)/sum(cos(as.numeric(entities_lats)*pi/180)))
  
  ## SIMULATION FULL 
  
  SPECTRA$SIM_full[[run]] <- list(
    ISOT = list(), 
    ITPC = list()
  )
  
  for(ii in 1:length(SPECTRA$entities_spec_rec)){
    print(ii)
    entity = SPECTRA$entities_spec_rec[ii]
    name = paste0("ENTITY", entity)
    data_yearly = DATA_past1000$CAVES$yearly_res[[run]] %>% filter(entity_id == entity)
    
    data  = ts(data = rev(data_yearly$ISOT), start = LastElement(data_yearly$year_BP), end = FirstElement(data_yearly$year_BP), deltat = 1)
    SPECTRA$SIM_full[[run]]$ISOT[[name]] = spectrum(data, plot = F)
    
    data  = ts(data = na.omit(rev(data_yearly$ITPC)), start = FirstElement(data_yearly$year_BP)-length(na.omit(data_yearly$ITPC)),
               end = FirstElement(data_yearly$year_BP), deltat = 1)
    SPECTRA$SIM_full[[run]]$ITPC[[name]] = spectrum(na.omit(data), plot = F)
  }
  
  SPECTRA$MEAN_SPEC$ISOT$full[[paste0("SIM_full_ISOT_", run)]] <- MeanSpectrum(SPECTRA$SIM_full[[run]]$ISOT)
  SPECTRA$MEAN_SPEC$ITPC$full[[paste0("SIM_full_ITPC_", run)]] <- MeanSpectrum(SPECTRA$SIM_full[[run]]$ITPC)
  
  SPECTRA$MEAN_SPEC_WEIGH$ISOT$full[[paste0("SIM_full_ISOT_", run)]] <- MeanSpectrum(SPECTRA$SIM_full[[run]]$ISOT, 
                                                                                     weights = cos(as.numeric(entities_lats)*pi/180)/sum(cos(as.numeric(entities_lats)*pi/180)))
  SPECTRA$MEAN_SPEC_WEIGH$ITPC$full[[paste0("SIM_full_ITPC_", run)]] <- MeanSpectrum(SPECTRA$SIM_full[[run]]$ITPC, 
                                                                                     weights = cos(as.numeric(entities_lats)*pi/180)/sum(cos(as.numeric(entities_lats)*pi/180)))
}

##Mean
for(var in c("ISOT","ITPC")){
  SPECTRA$MEAN_SPEC[[paste0("SIM_ds_", var)]] <- MeanSpectrum(list(SPECTRA$MEAN_SPEC[[var]]$ds[[paste0("SIM_ds_",var,"_a")]]$spec,
                                                                   SPECTRA$MEAN_SPEC[[var]]$ds[[paste0("SIM_ds_",var,"_b")]]$spec,
                                                                   SPECTRA$MEAN_SPEC[[var]]$ds[[paste0("SIM_ds_",var,"_c")]]$spec))
  SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_ds_", var)]] <- MeanSpectrum(list(SPECTRA$MEAN_SPEC_WEIGH[[var]]$ds[[paste0("SIM_ds_",var,"_a")]]$spec,
                                                                         SPECTRA$MEAN_SPEC_WEIGH[[var]]$ds[[paste0("SIM_ds_",var,"_b")]]$spec,
                                                                         SPECTRA$MEAN_SPEC_WEIGH[[var]]$ds[[paste0("SIM_ds_",var,"_c")]]$spec))
  SPECTRA$MEAN_SPEC[[paste0("SIM_full_", var)]] <- MeanSpectrum(list(SPECTRA$MEAN_SPEC[[var]]$full[[paste0("SIM_full_",var,"_a")]]$spec,
                                                                     SPECTRA$MEAN_SPEC[[var]]$full[[paste0("SIM_full_",var,"_b")]]$spec,
                                                                     SPECTRA$MEAN_SPEC[[var]]$full[[paste0("SIM_full_",var,"_c")]]$spec))
  SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_", var)]] <- MeanSpectrum(list(SPECTRA$MEAN_SPEC_WEIGH[[var]]$full[[paste0("SIM_full_",var,"_a")]]$spec,
                                                                           SPECTRA$MEAN_SPEC_WEIGH[[var]]$full[[paste0("SIM_full_",var,"_b")]]$spec,
                                                                           SPECTRA$MEAN_SPEC_WEIGH[[var]]$full[[paste0("SIM_full_",var,"_c")]]$spec))
}


#################################################
## Sim full FILTER Spectra ######################
# filtere full Sim signal auf down sim signal
# filtere full sim signal auf rec signal
# filtere down sim signal auf rec signal

source('Functions/Filter/EASY_Sensor_WM4.R')
source('Functions/Filter/filter_function3.R')
source('Functions/filter_window.R')


SPECTRA$SIM_filter_full_down <- list(ISOT = list(), ITPC = list())
SPECTRA$SIM_filter_full_rec <- list(ISOT = list(), ITPC = list())
SPECTRA$SIM_filter_down_rec <- list(ISOT = list(), ITPC = list())

filter <- list(
  ISOT = c(2.0, 16.0, 8.0),
  ITPC = c(3.0, 9.0, 3.0)
)


for(ii in 1:length(SPECTRA$entities_spec_rec)){
  print(ii)
  entity = SPECTRA$entities_spec_rec[ii]
  name = paste0("ENTITY", entity)
  
  
  for(run in c("a", "b", "c")){
    
    data_yearly = DATA_past1000$CAVES$yearly_res[[run]] %>% filter(entity_id == entity)
    
    #full -> down
    for(var in c("ISOT","ITPC")){
      # to perform easy_sensor_wm4 the data needs to go into the filter, forward in time... data_yearly comes in time from 1140BP -> 50BP so it's already in the right direction
      Results <- easy_sensor_wm4(1.0, na.omit(data_yearly[[var]]), filter[[var]][1])
      #the filtered Results are longer than the input, a new time has to be adjusted...
      time_new = seq(FirstElement(data_yearly$year_BP), LastElement(data_yearly$year_BP)-(length(Results)-length(data_yearly$ITPC)), by = -1)
      #since the first few years and the last include filter-artefacts, they are cut from the final output...
      # window nochmal neu einstellen!!!
      data  = filter_window(ts(data = rev(Results), start = LastElement(time_new), end = FirstElement(time_new), deltat = 1))
      SPECTRA$SIM_filter_full_down[[var]][[run]][[name]] = spectrum(data, plot = F)
    }
    
    #full-> rec
    for(var in c("ISOT","ITPC")){
      Results <- easy_sensor_wm4(1.0, na.omit(data_yearly[[var]]), filter[[var]][2])
      time_new = seq(FirstElement(data_yearly$year_BP), LastElement(data_yearly$year_BP)-(length(Results)-length(data_yearly$ITPC)), by = -1)
      data  = filter_window(ts(data = rev(Results), start = LastElement(time_new), end = FirstElement(time_new), deltat = 1))
      
      SPECTRA$SIM_filter_full_rec[[var]][[run]][[name]] = spectrum(data, plot = F)
    }
  }
  
  data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
  
  #down->rec

  diff = floor((LastElement(data_rec$interp_age)-FirstElement(data_rec$interp_age))/length(data_rec$interp_age))
  if(diff<1){diff = 1}
  
  for(run in c("a", "b", "c")){
    for(var in c("ISOT","ITPC")){
      record <- PaleoSpec::MakeEquidistant(data_rec$interp_age, data_rec[[paste0(var,"_", run)]],
                                           time.target = seq(from = FirstElement(data_rec$interp_age), to = LastElement(data_rec$interp_age), by = diff))
      # to perform easy_sensor_wm4 the data needs to go into the filter, forward in time... that's why the timeseries needs to be reversed
      Results <- easy_sensor_wm4(diff, na.omit(rev(record)), filter[[var]][3])
      time_new = seq(LastElement(data_rec$interp_age), LastElement(data_rec$interp_age)-diff*length(Results), by = -1*diff)
      data = filter_window(ts(data = rev(Results), start = LastElement(time_new), end = FirstElement(time_new), deltat = diff))
      SPECTRA$SIM_filter_down_rec[[var]][[run]][[name]] = spectrum(data, plot = F)
    }
  }
}

for(var in c("ISOT", "ITPC")){
  for(run in c("a", "b", "c")){
    SPECTRA$MEAN_SPEC[[var]][[paste0("full_down_", run)]] <- MeanSpectrum(SPECTRA$SIM_filter_full_down[[var]][[run]])
    SPECTRA$MEAN_SPEC_WEIGH[[var]][[paste0("full_down_", run)]] <- MeanSpectrum(SPECTRA$SIM_filter_full_down[[var]][[run]],
                                                                                weights = cos(as.numeric(entities_lats)*pi/180)/sum(cos(as.numeric(entities_lats)*pi/180)))
    SPECTRA$MEAN_SPEC[[var]][[paste0("full_rec_", run)]] <- MeanSpectrum(SPECTRA$SIM_filter_full_rec[[var]][[run]])
    SPECTRA$MEAN_SPEC_WEIGH[[var]][[paste0("full_rec_", run)]] <- MeanSpectrum(SPECTRA$SIM_filter_full_rec[[var]][[run]],
                                                                               weights = cos(as.numeric(entities_lats)*pi/180)/sum(cos(as.numeric(entities_lats)*pi/180)))
    SPECTRA$MEAN_SPEC[[var]][[paste0("down_rec_", run)]] <- MeanSpectrum(SPECTRA$SIM_filter_down_rec[[var]][[run]])
    SPECTRA$MEAN_SPEC_WEIGH[[var]][[paste0("down_rec_", run)]] <- MeanSpectrum(SPECTRA$SIM_filter_down_rec[[var]][[run]], 
                                                                               weights = cos(as.numeric(entities_lats)*pi/180)/sum(cos(as.numeric(entities_lats)*pi/180)))
  }
  
  SPECTRA$MEAN_SPEC[[paste0("SIM_full_down_",var)]]<- MeanSpectrum(list(SPECTRA$MEAN_SPEC[[var]]$full_down_a$spec, SPECTRA$MEAN_SPEC[[var]]$full_down_b$spec, SPECTRA$MEAN_SPEC[[var]]$full_down_c$spec))
  SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_down_",var)]]<- MeanSpectrum(list(SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_down_a$spec,SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_down_b$spec,SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_down_c$spec))
  SPECTRA$MEAN_SPEC[[paste0("SIM_full_rec_",var)]]<- MeanSpectrum(list(SPECTRA$MEAN_SPEC[[var]]$full_rec_a$spec, SPECTRA$MEAN_SPEC[[var]]$full_rec_b$spec,SPECTRA$MEAN_SPEC[[var]]$full_rec_c$spec))
  SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_rec_",var)]]<- MeanSpectrum(list(SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_rec_a$spec,SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_rec_b$spec,SPECTRA$MEAN_SPEC_WEIGH[[var]]$full_rec_c$spec))
  SPECTRA$MEAN_SPEC[[paste0("SIM_down_rec_",var)]]<- MeanSpectrum(list(SPECTRA$MEAN_SPEC[[var]]$down_rec_a$spec, SPECTRA$MEAN_SPEC[[var]]$down_rec_b$spec,SPECTRA$MEAN_SPEC[[var]]$down_rec_c$spec))
  SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_down_rec_",var)]]<- MeanSpectrum(list(SPECTRA$MEAN_SPEC_WEIGH[[var]]$down_rec_a$spec,SPECTRA$MEAN_SPEC_WEIGH[[var]]$down_rec_b$spec,SPECTRA$MEAN_SPEC_WEIGH[[var]]$down_rec_c$spec))
}

remove(entities_spec, diff, entity, ii, length, length_cave, name, record, Results, run, site, start_ts, stop_ts, var, easy_sensor_wm4, filter_function3,
       simpleawmean, simpleawsd, data_rec, data_yearly, entities_lats, filter, time_new, data)
