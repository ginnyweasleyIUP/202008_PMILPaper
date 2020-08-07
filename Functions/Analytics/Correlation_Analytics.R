library(dplyr)
library(latex2exp)
library(zoo)
library(PaleoSpec)
library(nest)
source("Functions/projection_ptlyr.R")
source("Functions/STACYmap_PMIL.R")
source("Functions/STACYmap_PMIL_NAgrid.R")

# this is correlation with downsampled temp and prec

#################################################
## ANALYSIS

No.digits = 2

CORR_CAVE <- list(entity_id = numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                  lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                  lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])))

for(run in c("a", "b", "c")){
  for(var in c("TEMP", "PREC", "ITPC")){
    CORR_CAVE[[paste0("CORR_",var, "_", run)]] <- numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var]))
    CORR_CAVE[[paste0("p_", var, "_", run)]] <- numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var]))
    CORR_CAVE[[paste0("ci1_", var, "_", run)]] <- numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var]))
    CORR_CAVE[[paste0("ci2_", var, "_", run)]] <- numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var]))
    for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])){
      entity = DATA_past1000$CAVES$entity_info$entity_id[mask_var][ii]
      CORR_CAVE$entity_id[ii] = entity
      CORR_CAVE$lon[ii] = DATA_past1000$CAVES$entity_info$longitude[mask_var][ii]
      CORR_CAVE$lat[ii] = DATA_past1000$CAVES$entity_info$latitude[mask_var][ii]
      data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
      # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
      double_time <- data_rec %>% group_by(interp_age) %>% count() %>% filter(n>1)
      data_rec <- data_rec %>% filter(!interp_age %in% double_time$interp_age)
      if(length(data_rec$interp_age)>4){
        record = zoo(x = data_rec$d18O_measurement, order.by = data_rec$interp_age)
        sim = zoo(x = data_rec[[paste0(var,"_", run)]], order.by = data_rec$interp_age)
        COR <- nexcf_ci(record, sim)
        
        CORR_CAVE[[paste0("CORR_",var, "_", run)]][ii] = COR$rxy
        CORR_CAVE[[paste0("p_",var, "_", run)]][ii] = COR$pval
        CORR_CAVE[[paste0("ci1_",var, "_", run)]][ii] = COR$ci.rxy[1]
        CORR_CAVE[[paste0("ci2_",var, "_", run)]][ii] = COR$ci.rxy[2]
        
      }else{
        CORR_CAVE[[paste0("CORR_",var, "_", run)]][ii] = NA
        CORR_CAVE[[paste0("p_",var, "_", run)]][ii] = NA
        CORR_CAVE[[paste0("ci1_",var, "_", run)]][ii] = NA
        CORR_CAVE[[paste0("ci2_",var, "_", run)]][ii] = NA
      }
    }
  }
}


print(paste0("Max Corr Temp: eID", CORR_CAVE$entity_id[CORR_CAVE$p_TEMP_a<0.1][which.max(abs(CORR_CAVE$CORR_TEMP_a[CORR_CAVE$p_TEMP_a<0.1]))], 
             " cave ", DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == CORR_CAVE$entity_id[CORR_CAVE$p_TEMP_a<0.1][which.max(abs(CORR_CAVE$CORR_TEMP_a[CORR_CAVE$p_TEMP_a<0.1]))]],
             ", rho = ", round(CORR_CAVE$CORR_TEMP_a[CORR_CAVE$p_TEMP_a<0.1][which.max(abs(CORR_CAVE$CORR_TEMP_a[CORR_CAVE$p_TEMP_a<0.1]))], digits = No.digits), 
             ", CI: (", round(CORR_CAVE$ci1_TEMP_a[CORR_CAVE$p_TEMP_a<0.1][which.max(abs(CORR_CAVE$CORR_TEMP_a[CORR_CAVE$p_TEMP_a<0.1]))], digits = No.digits),
             ", ", round(CORR_CAVE$ci2_TEMP_a[CORR_CAVE$p_TEMP_a<0.1][which.max(abs(CORR_CAVE$CORR_TEMP_a[CORR_CAVE$p_TEMP_a<0.1]))], digits = No.digits),")"))

print(paste0("Max Corr Precip: eID", CORR_CAVE$entity_id[CORR_CAVE$p_PREC_a<0.1][which.max(abs(CORR_CAVE$CORR_PREC_a[CORR_CAVE$p_PREC_a<0.1]))], 
             " cave ", DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == CORR_CAVE$entity_id[CORR_CAVE$p_PREC_a<0.1][which.max(abs(CORR_CAVE$CORR_PREC_a[CORR_CAVE$p_PREC_a<0.1]))]],
             ", rho = ", round(CORR_CAVE$CORR_PREC_a[CORR_CAVE$p_PREC_a<0.1][which.max(abs(CORR_CAVE$CORR_PREC_a[CORR_CAVE$p_PREC_a<0.1]))], digits = No.digits), 
             ", CI: (", round(CORR_CAVE$ci1_PREC_a[CORR_CAVE$p_PREC_a<0.1][which.max(abs(CORR_CAVE$CORR_PREC_a[CORR_CAVE$p_PREC_a<0.1]))], digits = No.digits),
             ", ", round(CORR_CAVE$ci2_PREC_a[CORR_CAVE$p_PREC_a<0.1][which.max(abs(CORR_CAVE$CORR_PREC_a[CORR_CAVE$p_PREC_a<0.1]))], digits = No.digits),")"))

print(paste0("Significant correlation for temperature: (", sum(CORR_CAVE$p_TEMP_a<0.1, na.rm = T), ", ",sum(CORR_CAVE$p_TEMP_b<0.1, na.rm = T), ", ", sum(CORR_CAVE$p_TEMP_c<0.1, na.rm = T), ")"))
print(paste0("Significant correlation for precipitation: (", sum(CORR_CAVE$p_PREC_a<0.1, na.rm = T), ", ",sum(CORR_CAVE$p_PREC_b<0.1, na.rm = T), ", ", sum(CORR_CAVE$p_PREC_c<0.1, na.rm = T), ")"))
print(paste0("Significant correlation for d18Opw: (", sum(CORR_CAVE$p_ITPC_a<0.1, na.rm = T), ", ",sum(CORR_CAVE$p_ITPC_b<0.1, na.rm = T), ", ", sum(CORR_CAVE$p_ITPC_c<0.1, na.rm = T), ")"))

rm(COR, CORR_CAVE, data_rec, double_time, entity, ii, No.digits, record, run, sim, var)
