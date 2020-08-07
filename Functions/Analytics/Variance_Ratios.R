# Variances for simulation and for proxy:
No.digits = 2

VARIANCE <- list(
  var_d18Oc = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_proxy_a = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_proxy_b = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_proxy_c = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_a = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_b = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_c = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_a = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_b = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_c = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_a_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_b_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_c_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_a_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_b_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_c_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
)

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
  data_yearly_a = DATA_past1000$CAVES$yearly_res$a %>% filter(entity_id == entity)
  data_yearly_b = DATA_past1000$CAVES$yearly_res$b %>% filter(entity_id == entity)
  data_yearly_c = DATA_past1000$CAVES$yearly_res$c %>% filter(entity_id == entity)
  VARIANCE$var_d18Oc[ii] = var(data_rec$d18O_measurement, na.rm = T)
  VARIANCE$var_proxy_a[ii] = var(data_rec$d18O_dw_eq_a, na.rm = T)
  VARIANCE$var_proxy_b[ii] = var(data_rec$d18O_dw_eq_b, na.rm = T)
  VARIANCE$var_proxy_c[ii] = var(data_rec$d18O_dw_eq_c, na.rm = T)
  VARIANCE$var_ds_a[ii] = var(data_rec$ITPC_a, na.rm = T)
  VARIANCE$var_ds_b[ii] = var(data_rec$ITPC_b, na.rm = T)
  VARIANCE$var_ds_c[ii] = var(data_rec$ITPC_c, na.rm = T)
  VARIANCE$var_full_a[ii] = var(data_yearly_a$ITPC, na.rm = T)
  VARIANCE$var_full_b[ii] = var(data_yearly_b$ITPC, na.rm = T)
  VARIANCE$var_full_c[ii] = var(data_yearly_c$ITPC, na.rm = T)
  VARIANCE$var_ds_a_isot[ii] = var(data_rec$ISOT_a, na.rm = T)
  VARIANCE$var_ds_b_isot[ii] = var(data_rec$ISOT_b, na.rm = T)
  VARIANCE$var_ds_c_isot[ii] = var(data_rec$ISOT_c, na.rm = T)
  VARIANCE$var_full_a_isot[ii] = var(data_yearly_a$ISOT, na.rm = T)
  VARIANCE$var_full_b_isot[ii] = var(data_yearly_b$ISOT, na.rm = T)
  VARIANCE$var_full_c_isot[ii] = var(data_yearly_c$ISOT, na.rm = T)
}

VARIANCE$var_proxy <- as.numeric(c(VARIANCE$var_proxy_a, VARIANCE$var_proxy_b, VARIANCE$var_proxy_c))
VARIANCE$var_ds <- as.numeric(c(VARIANCE$var_ds_a, VARIANCE$var_ds_b, VARIANCE$var_ds_c))
VARIANCE$var_full <- as.numeric(c(VARIANCE$var_full_a, VARIANCE$var_full_b, VARIANCE$var_full_c))
VARIANCE$var_ds_isot <- as.numeric(c(VARIANCE$var_ds_a_isot, VARIANCE$var_ds_b_isot, VARIANCE$var_ds_c_isot))
VARIANCE$var_full_isot <- as.numeric(c(VARIANCE$var_full_a_isot, VARIANCE$var_full_b_isot, VARIANCE$var_full_c_isot))


# calculate confidence interval via bootstrapping:

bstrap <- c()
for (i in 1:1000){
  bstrap <- c(bstrap,mean(sample((VARIANCE$var_d18Oc/VARIANCE$var_proxy_a)[mask_var],96,replace=T), na.rm = T))
}
print(paste0("VR between d18Oc and the drip water equivalent: ",
             round(mean((VARIANCE$var_d18Oc/VARIANCE$var_proxy_a)[mask_var], na.rm = T), digits = No.digits),", CI: (",
             round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.05), digits = No.digits),")"))


bstrap <- c()
for (i in 1:1000){
  bstrap <- c(bstrap,median(sample(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds")]][rep(mask_var,3)],96,replace=T), na.rm = T))
}
print(paste0("Median VR (down): ", round(median(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds")]][rep(mask_var,3)], na.rm = T), digits = No.digits),
             ", 90% CI: (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),")"))
bstrap <- c()
for (i in 1:1000){
  bstrap <- c(bstrap,median(sample(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full")]][rep(mask_var,3)],96,replace=T), na.rm = T))
}
print(paste0("Median VR (full): ", round(median(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full")]][rep(mask_var,3)], na.rm = T), digits = No.digits),
             ", 90% CI: (",round(quantile(bstrap,0.05), digits = No.digits),", ",round(quantile(bstrap,0.95), digits = No.digits),")"))


print(paste0("Max Outliner: eID", DATA_past1000$CAVES$entity_info$entity_id[mask_var][which.max((VARIANCE$var_proxy_a/VARIANCE$var_ds_a)[mask_var])], 
             " cave ", DATA_past1000$CAVES$entity_info$site_id[mask_var][which.max((VARIANCE$var_proxy_a/VARIANCE$var_ds_a)[mask_var])],
             ", VR = ", round((VARIANCE$var_proxy_a/VARIANCE$var_ds_a)[mask_var][which.max((VARIANCE$var_proxy_a/VARIANCE$var_ds_a)[mask_var])], digits = No.digits)))

print(paste0("Min Outliner: eID", DATA_past1000$CAVES$entity_info$entity_id[mask_var][which.min((VARIANCE$var_proxy_a/VARIANCE$var_ds_a)[mask_var])], 
             " cave ", DATA_past1000$CAVES$entity_info$site_id[mask_var][which.min((VARIANCE$var_proxy_a/VARIANCE$var_ds_a)[mask_var])],
             ", VR = ", round((VARIANCE$var_proxy_a/VARIANCE$var_ds_a)[mask_var][which.min((VARIANCE$var_proxy_a/VARIANCE$var_ds_a)[mask_var])], digits = No.digits)))

rm(data_rec, data_yearly_a, data_yearly_b, data_yearly_c, bstrap, entity)
rm(i, ii, No.digits)
