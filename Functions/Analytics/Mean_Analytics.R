# Mean State Analytics

library(plyr)
library(dplyr)
library(tidyverse)

No.digits = 2
load("Data/LM_HadCM3_annualmean.R")

for(run in c("a", "b", "c")){
  print(paste0("RUN ", run))
  
  mean_bias_full <- list()
  var_ratio_full <- list()
  mean_bias_ds <- list()
  var_ratio_ds <- list()
  mean_bias_lat <- list()
  for(entity in DATA_past1000$CAVES$entity_info$entity_id[mask_mean]){
    data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
    data_yearly = DATA_past1000$CAVES$yearly_res[[run]] %>% filter(entity_id == entity)
    
    mean_bias_full <- c(mean_bias_full, mean(data_yearly$ITPC, na.rm = T) - mean(data_rec[[paste0("d18O_dw_eq_", run)]], na.rm = T))
    mean_bias_ds <- c(mean_bias_ds, mean(data_rec[[paste0("ITPC_", run)]], na.rm = T) - mean(data_rec[[paste0("d18O_dw_eq_", run)]], na.rm = T))
    mean_bias_lat <- c(mean_bias_lat, DATA_past1000$CAVES$entity_info$latitude[DATA_past1000$CAVES$entity_info$entity_id == entity])
    var_ratio_full = c(var_ratio_full, var(data_rec[[paste0("d18O_dw_eq_", run)]], na.rm = T)/var(data_yearly$ITPC, na.rm = T))
    var_ratio_ds = c(var_ratio_ds, var(data_rec[[paste0("d18O_dw_eq_", run)]], na.rm = T)/var(data_rec[[paste0("ITPC_", run)]], na.rm = T))
  }
  
  weighing = cos(as.numeric(mean_bias_lat)*pi/180)/sum(cos(as.numeric(mean_bias_lat)*pi/180))
  
  if(method == "full"){
    c_left = mean(as.numeric(mean_bias_full), na.rm = T) - qnorm(0.95)*sd(as.numeric(mean_bias_full), na.rm = T)
    c_right = mean(as.numeric(mean_bias_full), na.rm = T) + qnorm(0.95)*sd(as.numeric(mean_bias_full), na.rm = T)
    print(paste0("Global Mean (full): ", round(mean(as.numeric(mean_bias_full), na.rm = T),digits = No.digits), "[‰],",
                 " 90% CI: (", round(c_left, digits = No.digits), ", ", round(c_right, digits = No.digits), ")"))  
    q_1 = (sum(weighing*as.numeric(mean_bias_full), na.rm = T)+ 1.644854*sqrt(sum(weighing*weighing)*sd(as.numeric(mean_bias_full), na.rm = T)))
    q_2 = (sum(weighing*as.numeric(mean_bias_full), na.rm = T)- 1.644854*sqrt(sum(weighing*weighing)*sd(as.numeric(mean_bias_full), na.rm = T)))
    print(paste0("Global Mean (full-weighted): ", round(sum(weighing*as.numeric(mean_bias_full), na.rm = T), digits = 3), " [‰],",
                 " 90% CI: (",round(q_2, digits = No.digits), ", ", round(q_1, digits = No.digits), ")"))
  } else{
    c_left = mean(as.numeric(mean_bias_ds), na.rm = T) - qnorm(0.95)*sd(as.numeric(mean_bias_ds), na.rm = T)
    c_right = mean(as.numeric(mean_bias_ds), na.rm = T) + qnorm(0.95)*sd(as.numeric(mean_bias_ds), na.rm = T)
    print(paste0("Global Mean (down): ", round(mean(as.numeric(mean_bias_ds), na.rm = T), digits = No.digits), "[‰],",
                 " 90% CI: (", round(c_left, digits = No.digits), ", ", round(c_right, digits = No.digits), ")"))
    
    q_1 = (sum(weighing*as.numeric(mean_bias_ds), na.rm = T)+ 1.644854*sqrt(sum(weighing*weighing)*sd(as.numeric(mean_bias_ds), na.rm = T)))
    q_2 = (sum(weighing*as.numeric(mean_bias_ds), na.rm = T)- 1.644854*sqrt(sum(weighing*weighing)*sd(as.numeric(mean_bias_ds), na.rm = T)))
    print(paste0("Global Mean (down-weighted): ", round(sum(weighing*as.numeric(mean_bias_ds), na.rm = T), digits = No.digits), " [‰],",
                 " 90% CI: (",round(q_2, digits = No.digits), ", ", round(q_1, digits = No.digits), ")"))
  }
  
  # Biggest Outliners:
  if(method == "full"){
    print(paste0("Max Outliner: eID", DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.max(as.numeric(mean_bias_full))], 
                 " cave ", DATA_past1000$CAVES$entity_info$site_id[mask_mean][which.max(as.numeric(mean_bias_full))],
                 ", Delta = ", round(as.numeric(mean_bias_full)[which.max(as.numeric(mean_bias_full))], digits = No.digits),
                 ", d18Oc_dw = ",round(mean(DATA_past1000$CAVES$record_res$d18O_dw_eq_a[DATA_past1000$CAVES$record_res$entity_id == DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.max(as.numeric(mean_bias_full))]],na.rm = T),digits = No.digits) , 
                 ", d18O = ", round(mean(DATA_past1000$CAVES$yearly_res$a$ITPC[DATA_past1000$CAVES$yearly_res$a$entity_id == DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.max(as.numeric(mean_bias_full))]],na.rm = T),digits = No.digits)))
    print(paste0("Min Outliner: eID", DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.min(as.numeric(mean_bias_full))], 
                 " cave ", DATA_past1000$CAVES$entity_info$site_id[mask_mean][which.min(as.numeric(mean_bias_full))],
                 ", Delta = ", round(as.numeric(mean_bias_full)[which.min(as.numeric(mean_bias_full))], digits = No.digits),
                 ", d18Oc_dw = ",round(mean(DATA_past1000$CAVES$record_res$d18O_dw_eq_a[DATA_past1000$CAVES$record_res$entity_id == DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.min(as.numeric(mean_bias_full))]],na.rm = T),digits = No.digits) , 
                 ", d18O = ", round(mean(DATA_past1000$CAVES$yearly_res$a$ITPC[DATA_past1000$CAVES$yearly_res$a$entity_id == DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.min(as.numeric(mean_bias_full))]],na.rm = T),digits = No.digits)))
  }else{
    print(paste0("Max Outliner: eID", DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.max(as.numeric(mean_bias_ds))], 
                 " cave ", DATA_past1000$CAVES$entity_info$site_id[mask_mean][which.max(as.numeric(mean_bias_full))],
                 ", Delta = ", round(as.numeric(mean_bias_ds)[which.max(as.numeric(mean_bias_ds))], digits = No.digits),
                 ", d18Oc_dw = ",round(mean(DATA_past1000$CAVES$record_res$d18O_dw_eq_a[DATA_past1000$CAVES$record_res$entity_id == DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.max(as.numeric(mean_bias_ds))]],na.rm = T),digits = No.digits) , 
                 ", d18O = ", round(mean(DATA_past1000$CAVES$record_res$ITPC_a[DATA_past1000$CAVES$record_res$entity_id == DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.max(as.numeric(mean_bias_ds))]],na.rm = T),digits = No.digits)))
    print(paste0("Min Outliner: eID", DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.min(as.numeric(mean_bias_ds))], 
                 " cave ", DATA_past1000$CAVES$entity_info$site_id[mask_mean][which.min(as.numeric(mean_bias_full))],
                 ", Delta = ", round(as.numeric(mean_bias_ds)[which.min(as.numeric(mean_bias_ds))], digits = No.digits),
                 ", d18Oc_dw = ",round(mean(DATA_past1000$CAVES$record_res$d18O_dw_eq_a[DATA_past1000$CAVES$record_res$entity_id == DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.min(as.numeric(mean_bias_ds))]],na.rm = T),digits = No.digits) , 
                 ", d18O = ", round(mean(DATA_past1000$CAVES$record_res$ITPC_a[DATA_past1000$CAVES$record_res$entity_id == DATA_past1000$CAVES$entity_info$entity_id[mask_mean][which.min(as.numeric(mean_bias_ds))]],na.rm = T),digits = No.digits)))
  }
  
  if(method == "full"){
    COR <- cor.test(as.numeric(mean_bias_full)[mask_var[mask_mean]], as.numeric(var_ratio_full)[mask_var[mask_mean]], conf.level = 0.9)
    print(paste0("Correlation between mean bias and var ratio (full): r=", round(COR$estimate, digits = No.digits), 
                 " CI: (", round(COR$conf.int[1], digits = No.digits), ", ",round(COR$conf.int[2], digits = No.digits), "),",
                 " p-level: ", round(COR$p.value, digits = No.digits)))
  }else{
    COR <- cor.test(as.numeric(mean_bias_ds)[mask_var[mask_mean]], as.numeric(var_ratio_ds)[mask_var[mask_mean]], conf.level = 0.9)
    print(paste0("Correlation between mean bias and var ratio (ds): r=", round(COR$estimate, digits = No.digits), 
                 " CI: (", round(COR$conf.int[1], digits = No.digits), ", ",round(COR$conf.int[2], digits = No.digits), "),",
                 " p-level: ", round(COR$p.value, digits = No.digits)))
    
  }
  
  # CLUSTER
  
  cluster_mean <- list()
  
  for(cluster in 1:9){
    entity_list <- as.tibble(DATA_past1000$CAVES$cluster_list) %>% filter(cluster_id == cluster)
    entity_list <- entity_list$entity_id
    
    mean_bias_full <- list()
    mean_bias_ds <- list()
    for(entity in entity_list){
      data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
      data_yearly = DATA_past1000$CAVES$yearly_res[[run]] %>% filter(entity_id == entity)
      mean_bias_full <- c(mean_bias_full, mean(data_yearly$ITPC, na.rm = T) - mean(data_rec[[paste0("d18O_dw_eq_", run)]], na.rm = T))
      mean_bias_ds <- c(mean_bias_ds,     mean(data_rec[[paste0("ITPC_", run)]], na.rm = T) - mean(data_rec[[paste0("d18O_dw_eq_", run)]], na.rm = T))
    }
    if(method == "full"){
      c_left = mean(as.numeric(mean_bias_full), na.rm = T) - qnorm(0.95)*sd(as.numeric(mean_bias_full), na.rm = T)
      c_right = mean(as.numeric(mean_bias_full), na.rm = T) + qnorm(0.95)*sd(as.numeric(mean_bias_full), na.rm = T)
      print(paste0("Cluster ",cluster," Mean (full): ", round(mean(as.numeric(mean_bias_full), na.rm = T),digits = No.digits), "[‰],",
                   " 90% CI: (", round(c_left, digits = No.digits), ", ", round(c_right, digits = No.digits), ")"))  
    }else{
      c_left = mean(as.numeric(mean_bias_ds), na.rm = T) - qnorm(0.95)*sd(as.numeric(mean_bias_ds), na.rm = T)
      c_right = mean(as.numeric(mean_bias_ds), na.rm = T) + qnorm(0.95)*sd(as.numeric(mean_bias_ds), na.rm = T)
      print(paste0("Cluster ",cluster," Mean (down): ", round(mean(as.numeric(mean_bias_ds), na.rm = T), digits = No.digits), "[‰], ",
                   "90% CI: (", round(c_left, digits = No.digits), ", ", round(c_right, digits = No.digits), ")"))
    }
    
    
  }
  
  
  
  sd.global <- apply(DATA_past1000$SIM_yearly_a$ITPC, c(1,2), sd, na.rm = T)
  weighing = array(dim = c(96,73))
  for(ii in 1:96){
    lats = seq(90,-90,length.out = 73)
    weighing[ii,] = cos(lats*pi/180)/sum(cos(lats*pi/180))/96
  }# cos(lats*pi/180)/sum(cos(lats*pi/180))
  q_1 = (sum(weighing*as.numeric(sd.global), na.rm = T)+ 1.644854*sqrt(sum(weighing*weighing)*sd(as.numeric(sd.global), na.rm = T)))
  q_2 = (sum(weighing*as.numeric(sd.global), na.rm = T)- 1.644854*sqrt(sum(weighing*weighing)*sd(as.numeric(sd.global), na.rm = T)))
  print(paste0("sd ITPC of simulation (full-weighted): ", round(sum(as.numeric(weighing*sd.global), na.rm = T), digits = 3), " [‰],",
               " 90% CI: (",round(q_2, digits = No.digits), ", ", round(q_1, digits = No.digits), ")"))
  
  # Temperature range for caves:
  
  print(paste0("simulated temperature range for caves: (", round(range(DATA_past1000$CAVES$entity_info$temp_ds)[1], digits = 3),
               ", ",round(range(DATA_past1000$CAVES$entity_info$temp_ds)[2], digits = 3) ,") [°C]"))
  
}












rm(mean_bias_ds, mean_bias_full, mean_bias_lat, entity_list, cluster, entity, data_rec, data_yearly_a, data_yearly_b, data_yearly_c, q_1, q_2, weighing, cluster_mean)
rm(var_ratio_ds, var_ratio_full, COR, No.digits, sd.global)
