#################################################
## Ensemble Analytics

source("Functions/network_mask.R")
source("Functions/network_mask_ensemble.R")

##make table in supplemente with exact numbers for correlations and SNR,
# only have percentage of increase in text and refer to supplements

No.digits = 2

#On a site level, we increased the tuning from a median of correlations of $c = 0.15 (0.03, 0.47)$ to $c = 0.51 (0.41,0.91)$.
list <- DATA_past1000$CAVES$entity_info %>% select(site_id, entity_id) %>%
  filter(entity_id %in% DATA_past1000$CAVES$entity_info$entity_id[mask_spec]) %>% group_by(site_id) %>% count() %>% filter(n>1)
CORR = NETWORK$record
CORR[NETWORK$record_p>0.1] = NA
CORR = CORR[network_mask(list,"site")]
CORR_tuned = c_ensemble[network_mask_ensemble(list, "site")]
SNR = abs(median(na.omit(as.numeric(CORR)), na.rm = T))/(1-abs(median(na.omit(as.numeric(CORR)), na.rm = T)))
SNR_tuned = abs(median(na.omit(as.numeric(CORR_tuned)), na.rm = T))/(1-abs(median(na.omit(as.numeric(CORR_tuned)), na.rm = T)))

bstrap <- c()
for (i in 1:1000){
  bstrap <- c(bstrap,median(sample(na.omit(as.numeric(CORR)),length(na.omit(as.numeric(CORR))),replace=T), na.rm = T))
}
bstrap_tuned <- c()
for (i in 1:1000){
  bstrap_tuned <- c(bstrap_tuned,median(sample(na.omit(as.numeric(CORR_tuned)),length(na.omit(as.numeric(CORR_tuned))),replace=T), na.rm = T))
}
print(paste0("On site level: median of correlations of c = ", round(median(na.omit(as.numeric(CORR)), na.rm = T), digits = No.digits),
             " (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),
             "), after tuning: c = ",round(median(na.omit(as.numeric(CORR_tuned)), na.rm = T), digits = No.digits),
             " (",round(quantile(bstrap_tuned,0.05), digits = No.digits),", ",round(quantile(bstrap_tuned,0.95), digits = No.digits),
             ")"))
print(paste0("Increase by ", round((median(na.omit(as.numeric(CORR_tuned)), na.rm = T)/median(na.omit(as.numeric(CORR)), na.rm = T)-1)*100, digits = No.digits),"%"))

bstrap <- c()
for (i in 1:1000){
  bsample = sample(na.omit(as.numeric(CORR)),length(na.omit(as.numeric(CORR))),replace=T)
  bstrap <- c(bstrap,abs(median(bsample, na.rm = T))/(1-abs(median(bsample, na.rm = T))))
}
bstrap_tuned <- c()
for (i in 1:1000){
  bsample = sample(na.omit(as.numeric(CORR_tuned)),length(na.omit(as.numeric(CORR_tuned))),replace=T)
  bstrap_tuned <- c(bstrap_tuned,abs(median(bsample, na.rm = T))/(1-abs(median(bsample, na.rm = T))))
}
print(paste0("On site level: SNR = ", round(SNR, digits = No.digits),
             " (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),
             "), after tuning: SNR = ",round(SNR_tuned, digits = No.digits),
             " (",round(quantile(bstrap_tuned,0.05), digits = No.digits),", ",round(quantile(bstrap_tuned,0.95), digits = No.digits),
             ")"))
print(paste0("Increase by ", round((SNR_tuned/SNR-1)*100, digits = No.digits),"%"))



#On a gridbox level, we obtained stronger but not significant correlation with changes from$c=0.14(-0.33,0.36)$ to $c=0.43 (-0.92, 0.92)$

list <- as.tibble(DATA_past1000$CAVES$gridbox_list) %>% group_by(gridbox_id) %>% count() %>% filter(n>1)
CORR = NETWORK$record
CORR[NETWORK$record_p>0.1] = NA
CORR = CORR[network_mask(list, "gridbox")]
CORR_tuned = c_ensemble[network_mask_ensemble(list, "gridbox")]
SNR = abs(median(na.omit(as.numeric(CORR)), na.rm = T))/(1-abs(median(na.omit(as.numeric(CORR)), na.rm = T)))
SNR_tuned = abs(median(na.omit(as.numeric(CORR_tuned)), na.rm = T))/(1-abs(median(na.omit(as.numeric(CORR_tuned)), na.rm = T)))

bstrap <- c()
for (i in 1:1000){
  bstrap <- c(bstrap,median(sample(na.omit(as.numeric(CORR)),length(na.omit(as.numeric(CORR))),replace=T), na.rm = T))
}
bstrap_tuned <- c()
for (i in 1:1000){
  bstrap_tuned <- c(bstrap_tuned,median(sample(na.omit(as.numeric(CORR_tuned)),length(na.omit(as.numeric(CORR_tuned))),replace=T), na.rm = T))
}
print(paste0("On gridbox level: median of correlations of c = ", round(median(na.omit(as.numeric(CORR)), na.rm = T), digits = No.digits),
             " (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),
             "), after tuning: c = ",round(median(na.omit(as.numeric(CORR_tuned)), na.rm = T), digits = No.digits),
             " (",round(quantile(bstrap_tuned,0.05), digits = No.digits),", ",round(quantile(bstrap_tuned,0.95), digits = No.digits),
             ")"))
print(paste0("Increase by ", round((median(na.omit(as.numeric(CORR_tuned)), na.rm = T)/median(na.omit(as.numeric(CORR)), na.rm = T)-1)*100, digits = No.digits),"%"))

bstrap <- c()
for (i in 1:1000){
  bsample = sample(na.omit(as.numeric(CORR)),length(na.omit(as.numeric(CORR))),replace=T)
  bstrap <- c(bstrap,abs(median(bsample, na.rm = T))/(1-abs(median(bsample, na.rm = T))))
}
bstrap_tuned <- c()
for (i in 1:1000){
  bsample = sample(na.omit(as.numeric(CORR_tuned)),length(na.omit(as.numeric(CORR_tuned))),replace=T)
  bstrap_tuned <- c(bstrap_tuned,abs(median(bsample, na.rm = T))/(1-abs(median(bsample, na.rm = T))))
}
print(paste0("On gridbox level: SNR = ", round(SNR, digits = No.digits),
             " (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),
             "), after tuning: SNR = ",round(SNR_tuned, digits = No.digits),
             " (",round(quantile(bstrap_tuned,0.05), digits = No.digits),", ",round(quantile(bstrap_tuned,0.95), digits = No.digits),
             ")"))
print(paste0("Increase by ", round((SNR_tuned/SNR-1)*100, digits = No.digits),"%"))

#################################################
# CLUSTER

cluster_list <- as.tibble(DATA_past1000$CAVES$cluster_list) %>% group_by(cluster_id) %>% count() %>% filter(n>1)
for(cluster in 1:9){
  
  list = as.tibble(DATA_past1000$CAVES$cluster_list) %>% filter(cluster_id == cluster)
  CORR = NETWORK$record
  CORR[NETWORK$record_p>0.1 | is.na(NETWORK$record_gauss_p)] = NA
  CORR = CORR[network_mask(list, "cluster")]
  CORR_tuned = c_ensemble[network_mask_ensemble(list, "cluster")]
  
  SNR = abs(median(na.omit(as.numeric(CORR)), na.rm = T))/(1-abs(median(na.omit(as.numeric(CORR)), na.rm = T)))
  SNR_tuned = abs(median(na.omit(as.numeric(CORR_tuned)), na.rm = T))/(1-abs(median(na.omit(as.numeric(CORR_tuned)), na.rm = T)))
  
  if(all(is.na(CORR))){
    print(paste0("c", cluster, " all NA"))
  }else{
    bstrap <- c()
    for (i in 1:1000){
      bstrap <- c(bstrap,median(sample(na.omit(as.numeric(CORR)),length(na.omit(as.numeric(CORR))),replace=T), na.rm = T))
    }
    print(paste0("On cluster",cluster," level: median of correlations of c = ", round(median(na.omit(as.numeric(CORR)), na.rm = T), digits = No.digits),
                 " (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),
                 ")"))
    bstrap <- c()
    for (i in 1:1000){
      bsample = sample(na.omit(as.numeric(CORR)),length(na.omit(as.numeric(CORR))),replace=T)
      bstrap <- c(bstrap,abs(median(bsample, na.rm = T))/(1-abs(median(bsample, na.rm = T))))
    }
    print(paste0("On cluster",cluster," level: SNR = ", round(SNR, digits = No.digits),
                 " (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),
                 ")"))
  }

  if(all(is.na(CORR_tuned))){
    print(paste0("c_tuned", cluster, " all NA"))
  }else{
    bstrap_tuned <- c()
    for (i in 1:1000){
      bstrap_tuned <- c(bstrap_tuned,median(sample(na.omit(as.numeric(CORR_tuned)),length(na.omit(as.numeric(CORR_tuned))),replace=T), na.rm = T))
    }
    print(paste0("   after tuning: c = ",round(median(na.omit(as.numeric(CORR_tuned)), na.rm = T), digits = No.digits),
                 " (",round(quantile(bstrap_tuned,0.05), digits = No.digits),", ",round(quantile(bstrap_tuned,0.95), digits = No.digits),
                 ")"))
    
    bstrap_tuned <- c()
    for (i in 1:1000){
      bsample = sample(na.omit(as.numeric(CORR_tuned)),length(na.omit(as.numeric(CORR_tuned))),replace=T)
      bstrap_tuned <- c(bstrap_tuned,abs(median(bsample, na.rm = T))/(1-abs(median(bsample, na.rm = T))))
    }
    print(paste0("   after tuning: SNR = ",round(SNR_tuned, digits = No.digits),
                 " (",round(quantile(bstrap_tuned,0.05), digits = No.digits),", ",round(quantile(bstrap_tuned,0.95), digits = No.digits),
                 ")"))
  }
  
  
}

#################################################

# GLOBAL
CORR = NETWORK$record
CORR[NETWORK$record_p>0.1 | is.na(NETWORK$record_gauss_p)] = NA
CORR_tuned = c_ensemble
SNR = abs(median(na.omit(as.numeric(CORR)), na.rm = T))/(1-abs(median(na.omit(as.numeric(CORR)), na.rm = T)))
SNR_tuned = abs(median(na.omit(as.numeric(CORR_tuned)), na.rm = T))/(1-abs(median(na.omit(as.numeric(CORR_tuned)), na.rm = T)))

bstrap <- c()
for (i in 1:1000){
  bstrap <- c(bstrap,median(sample(na.omit(as.numeric(CORR)),length(na.omit(as.numeric(CORR))),replace=T), na.rm = T))
}
bstrap_tuned <- c()
for (i in 1:1000){
  bstrap_tuned <- c(bstrap_tuned,median(sample(na.omit(as.numeric(CORR_tuned)),length(na.omit(as.numeric(CORR_tuned))),replace=T), na.rm = T))
}
print(paste0("On global level: median of correlations of c = ", round(median(na.omit(as.numeric(CORR)), na.rm = T), digits = No.digits),
             " (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),
             "), after tuning: c = ",round(median(na.omit(as.numeric(CORR_tuned)), na.rm = T), digits = No.digits),
             " (",round(quantile(bstrap_tuned,0.05), digits = No.digits),", ",round(quantile(bstrap_tuned,0.95), digits = No.digits),
             ")"))
print(paste0("Increase by ", round((median(na.omit(as.numeric(CORR_tuned)), na.rm = T)/median(na.omit(as.numeric(CORR)), na.rm = T)-1)*100, digits = No.digits),"%"))

bstrap <- c()
for (i in 1:1000){
  bsample = sample(na.omit(as.numeric(CORR)),length(na.omit(as.numeric(CORR))),replace=T)
  bstrap <- c(bstrap,abs(median(bsample, na.rm = T))/(1-abs(median(bsample, na.rm = T))))
}
bstrap_tuned <- c()
for (i in 1:1000){
  bsample = sample(na.omit(as.numeric(CORR_tuned)),length(na.omit(as.numeric(CORR_tuned))),replace=T)
  bstrap_tuned <- c(bstrap_tuned,abs(median(bsample, na.rm = T))/(1-abs(median(bsample, na.rm = T))))
}
print(paste0("On global level: SNR = ", round(SNR, digits = No.digits),
             " (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),
             "), after tuning: SNR = ",round(SNR_tuned, digits = No.digits),
             " (",round(quantile(bstrap_tuned,0.05), digits = No.digits),", ",round(quantile(bstrap_tuned,0.95), digits = No.digits),
             ")"))
print(paste0("Increase by ", round((SNR_tuned/SNR-1)*100, digits = No.digits),"%"))


rm(No.digits, bstrap, bstrap_tuned, CORR, CORR_tuned, i, SNR_tuned, SNR, list, cluster_list, bsample, cluster)
