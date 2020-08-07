#################################################
## Paper Figure 8 ###############################
#################################################

## Network 

#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(nest)
library(PaleoSpec)


lats = DATA_past1000$CAVES$entity_info$latitude[mask_spec]
longs = DATA_past1000$CAVES$entity_info$longitude[mask_spec]

dist<-fossil::earth.dist(cbind(lats,longs),dist=TRUE)
dist_matrix <- as.matrix(dist)

load("C_ensemble_max.RData")
c_ensemble <- C
ensemble_entities <- c(14,21,33,48,51,76,85,93,94,95,97,111,113,117,118,123,128,137,147,165,172,178,179,185,187,
                       201,202,209,212,222,226,238,240,242,253,278,286,289,294,296,305,326,329,330,335,351,361,378,390,395,399,
                       420,422,430,435,437,442,443,447,448,461,464,496,498,506,514,522,528,538,539,541,544,546,547,548,560,564,573,577,588,589,595,598,
                       613,620,621,623,672,673)
lat_e = numeric(length(ensemble_entities))
long_e = numeric(length(ensemble_entities))
for(ii in 1:length(ensemble_entities)){
  entity = ensemble_entities[ii]
  if(entity %in% DATA_past1000$CAVES$entity_info$entity_id){
    lat_e[ii] = DATA_past1000$CAVES$entity_info$latitude[DATA_past1000$CAVES$entity_info$entity_id == entity]
    long_e[ii] = DATA_past1000$CAVES$entity_info$longitude[DATA_past1000$CAVES$entity_info$entity_id == entity]
  }else{
    lat_e[ii] = NA
    long_e[ii] = NA
  }

}
dist_e<-fossil::earth.dist(cbind(lat_e,long_e),dist=TRUE)
dist_matrix_e <- as.matrix(dist_e)
rm(C)

#################################################

# NETWORK MATRIX

length_eID = length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])
NETWORK <- list(record = array(dim = c(length_eID, length_eID)), record_p = array(dim = c(length_eID, length_eID)),
                record_gauss = array(dim = c(length_eID, length_eID)), record_gauss_p = array(dim = c(length_eID, length_eID)),
                sim_full_a = array(dim = c(length_eID, length_eID)),       sim_full_a_p = array(dim = c(length_eID, length_eID)),
                sim_down_a = array(dim = c(length_eID, length_eID)),       sim_down_a_p = array(dim = c(length_eID, length_eID)),
                sim_full_gauss_a = array(dim = c(length_eID, length_eID)), sim_full_gauss_a_p = array(dim = c(length_eID, length_eID)),
                sim_down_gauss_a = array(dim = c(length_eID, length_eID)), sim_down_gauss_a_p = array(dim = c(length_eID, length_eID)),
                sim_full_b = array(dim = c(length_eID, length_eID)),       sim_full_b_p = array(dim = c(length_eID, length_eID)),
                sim_down_b = array(dim = c(length_eID, length_eID)),       sim_down_b_p = array(dim = c(length_eID, length_eID)),
                sim_full_gauss_b = array(dim = c(length_eID, length_eID)), sim_full_gauss_b_p = array(dim = c(length_eID, length_eID)),
                sim_down_gauss_b = array(dim = c(length_eID, length_eID)), sim_down_gauss_b_p = array(dim = c(length_eID, length_eID)),
                sim_full_c = array(dim = c(length_eID, length_eID)),       sim_full_c_p = array(dim = c(length_eID, length_eID)),
                sim_down_c = array(dim = c(length_eID, length_eID)),       sim_down_c_p = array(dim = c(length_eID, length_eID)),
                sim_full_gauss_c = array(dim = c(length_eID, length_eID)), sim_full_gauss_c_p = array(dim = c(length_eID, length_eID)),
                sim_down_gauss_c = array(dim = c(length_eID, length_eID)), sim_down_gauss_c_p = array(dim = c(length_eID, length_eID)))


print("record")
TS = list()
TS_gauss = list()

for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_spec]){
  data = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity) %>% select(interp_age, d18O_measurement) %>% rename(time = interp_age, value = d18O_measurement)

  # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
  double_time = data %>% group_by(time) %>% count() %>% filter(n>1)
  data = data %>% filter(!time %in% double_time$time) %>% filter(!is.na(value))
  TS[[paste0("Entity", entity)]] <- zoo(x = data$value, order.by = data$time)
  TS_gauss[[paste0("Entity", entity)]] <- gaussdetr(zoo(x = data$value, order.by = data$time), tsc.in = 100)$Xsmooth
}

C<-matrix(NA,nrow=length(TS),ncol=length(TS))
colnames(C)<-rownames(C)<-names(TS)
C_gauss <- P_gauss <- P <- C

for (i in 1:(length(TS)-1)){
  for (j in (i+1):length(TS)){
    temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
    temp_gauss<-nest::nexcf_ci(TS_gauss[[i]],TS_gauss[[j]],conflevel=0.1)
    C[i,j]<-temp$rxy
    P[i,j]<-P[j,i]<-temp$pval
    C[j,i]=C[i,j]
    rm(temp)
    C_gauss[i,j]<-temp_gauss$rxy
    P_gauss[i,j]<-P_gauss[j,i]<-temp_gauss$pval
    C_gauss[j,i]=C_gauss[i,j]
    rm(temp_gauss)
  }
}

NETWORK$record = C
NETWORK$record_p = P
NETWORK$record_gauss = C_gauss
NETWORK$record_gauss_p = P_gauss


for(run in c("a", "b", "c")){
#for(run in c("a")){  
  print(run)
  for(method in c("sim_full", "sim_down")){
    print(method)

    TS = list()
    TS_gauss = list()

    for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_spec]){
      if(method == "sim_down"){
        data = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity) %>% select(interp_age, ITPC_a) %>% rename(time = interp_age, value = ITPC_a)
      }else{
        data = DATA_past1000$CAVES$yearly_res$a %>% filter(entity_id == entity) %>% select(year_BP, ITPC) %>% rename (time = year_BP, value = ITPC)
      }

      # zoo cannot handle objects where order.by has two elements which is why they are sorted out here (no better option found)
      double_time = data %>% group_by(time) %>% count() %>% filter(n>1)
      data = data %>% filter(!time %in% double_time$time) %>% filter(!is.na(value))
      TS[[paste0("Entity", entity)]] <- zoo(x = data$value, order.by = data$time)
      TS_gauss[[paste0("Entity", entity)]] <- gaussdetr(zoo(x = data$value, order.by = data$time), tsc.in = 100)$Xsmooth
    }

    C<-matrix(NA,nrow=length(TS),ncol=length(TS))
    colnames(C)<-rownames(C)<-names(TS)
    C_gauss <- P_gauss <- P <- C

    for (i in 1:(length(TS)-1)){
      for (j in (i+1):length(TS)){
        temp<-nest::nexcf_ci(TS[[i]],TS[[j]],conflevel=0.1)
        temp_gauss<-nest::nexcf_ci(TS_gauss[[i]],TS_gauss[[j]],conflevel=0.1)
        C[i,j]<-temp$rxy
        P[i,j]<-P[j,i]<-temp$pval
        C[j,i]=C[i,j]
        rm(temp)
        C_gauss[i,j]<-temp_gauss$rxy
        P_gauss[i,j]<-P_gauss[j,i]<-temp_gauss$pval
        C_gauss[j,i]=C_gauss[i,j]
        rm(temp_gauss)
      }
    }

    NETWORK[[paste0(method,"_",run)]] = C
    NETWORK[[paste0(method,"_",run,"_p")]] = P
    NETWORK[[paste0(method,"_",run, "_gauss")]] = C_gauss
    NETWORK[[paste0(method,"_",run, "_gauss_p")]] = P_gauss

  }
}


#################################################

No.digits = 2
link_density = 0.05

C_SIM = NETWORK$sim_down_a
C_REC = NETWORK$record

o_sim = order(abs(C_SIM), na.last = F)
o_rec = order(abs(C_REC), na.last = F)
C_SIM[o_sim[1:floor((length(o_sim)-link_density*length(o_sim)))]] =  NA
C_REC[o_rec[1:floor((length(o_rec)-link_density*length(o_rec)))]] = NA

C_SIM[NETWORK$sim_down_a_p>0.1] = NA
C_REC[NETWORK$record_p>0.1] = NA

bstrap <- c()
for (i in 1:1000){
  bstrap <- c(bstrap,mean(sample(abs(na.omit(as.numeric(C_REC))),length(na.omit(as.numeric(C_REC))),replace=T), na.rm = T))
}
print(paste0("Absolut Mean 5% strongest (record): ", round(mean(abs(na.omit(as.numeric(C_REC))), na.rm = T), digits = No.digits),
             ", 90% CI: (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),")"))


bstrap <- c()
for (i in 1:1000){
  bstrap <- c(bstrap,mean(sample(abs(na.omit(as.numeric(C_SIM))),length(na.omit(as.numeric(C_SIM))),replace=T), na.rm = T))
}
print(paste0("Absolut Mean 5% strongest (sim): ", round(mean(abs(na.omit(as.numeric(C_SIM))), na.rm = T), digits = No.digits),
             ", 90% CI: (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),")"))

rm(C, C_gauss, data, P, P_gauss, TS, TS_gauss, dist, dist_e, ensemble_entities, entity, i, ii, j, length_eID, method, run)
rm(C_REC, C_SIM, bstrap, lat_e, lats, link_density, long_e, longs, No.digits, o_rec, o_sim)
