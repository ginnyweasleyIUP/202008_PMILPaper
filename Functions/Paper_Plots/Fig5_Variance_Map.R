#################################################
## Paper Figure 1 ###############################
#################################################

#################################################

# Variances for simulation and for proxy:

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
  var_full_c = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
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
}
  

## PLOTTING #####################################

library(latex2exp)

## Variance Map #################################

source("Functions/var_map_plot.R")

for(run in c("a","b","c")){
  Point_Lyr <- list(
    lon = DATA_past1000$CAVES$entity_info$longitude[mask_var],
    lat = DATA_past1000$CAVES$entity_info$latitude[mask_var],
    value = log10(VARIANCE$var_proxy_a[mask_var]/VARIANCE$var_ds_a[mask_var])
  )
  
  Point_Lyr$lon = Point_Lyr$lon[!is.na(Point_Lyr$value)]
  Point_Lyr$lat = Point_Lyr$lat[!is.na(Point_Lyr$value)]
  Point_Lyr$value = Point_Lyr$value[!is.na(Point_Lyr$value)]
  
  GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 2.5
  
  # Points in China werden gejittered und sind nicht an ihrer richtigen Position!
  plot_var <- var_map_plot(Point_Lyr =  Point_Lyr, pt_size =  3, txt_size =  10)
  plot_var
  
  plot_var %>% ggsave(filename = paste0('Fig5_Variance_1_map_xnap',run, '.pdf'), plot = ., path = 'Paper_Plots', 
                      width = 2*6, height = 2*PLOTTING_VARIABLES$HEIGHT/1.5, units = 'cm', dpi = 'print', device = "pdf")
  plot_var %>% ggsave(filename = paste0('Fig5_Variance_1_map_xnap',run, '.png'), plot = ., path = 'Paper_Plots', 
                      width = 2*6, height = 2*PLOTTING_VARIABLES$HEIGHT/1.5, units = 'cm', dpi = 'print', device = "png")
}

rm(VARIANCE, Point_Lyr, plot_var, data_rec, data_yearly_a, data_yearly_b, data_yearly_c, entity, ii, run)
