library(dplyr)
library(latex2exp)
library(zoo)
library(nest)
library(PaleoSpec)
source("Functions/projection_ptlyr.R")
source("Functions/STACYmap_PMIL.R")
source("Functions/STACYmap_PMIL_NAgrid.R")

# this is correlation with downsampled temp and prec

#################################################
## ANALYSIS

run = "a"
load("Data/LM_HadCM3_correlation.R")

CORR_FIELD <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)),
                   TEMP_P = array(dim = c(96,73)), PREC_P = array(dim = c(96,73)))

CORR_FIELD$TEMP = DATA_EXPORT_CORR[[run]]$TEMP
CORR_FIELD$TEMP_p = DATA_EXPORT_CORR[[run]]$TEMP_P
CORR_FIELD$PREC = DATA_EXPORT_CORR[[run]]$PREC
CORR_FIELD$PREC_p = DATA_EXPORT_CORR[[run]]$PREC_P

rm(DATA_EXPORT_CORR)

CORR_CAVE <- list(entity_id = numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                  lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                  lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])))

for(run in c("a")){
  for(var in c("TEMP", "PREC")){
    CORR_CAVE[[paste0("CORR_",var)]] <- numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var]))
    CORR_CAVE[[paste0("p_", var)]] <- numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var]))
    for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])){
      print(ii)
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
        
        CORR_CAVE[[paste0("CORR_",var)]][ii] = COR$rxy
        CORR_CAVE[[paste0("p_",var)]][ii] = COR$pval
      }else{
        CORR_CAVE[[paste0("CORR_",var)]][ii] = NA
        CORR_CAVE[[paste0("p_",var)]][ii] = NA
      }
    }
  }
}


#################################################
## PLOTS ########################################
#################################################
for(run in c("a")){
  print(run)
  Plot_lyr_temp = CORR_FIELD$TEMP
  Plot_lyr_temp_p =CORR_FIELD$TEMP_P
  Plot_lyr_temp[Plot_lyr_temp_p > 0.1] = NA
  Plot_lyr_prec = CORR_FIELD$PREC
  Plot_lyr_prec_p = CORR_FIELD$PREC_P
  Plot_lyr_prec[Plot_lyr_prec_p > 0.1] = NA
  
  Plot_lyr_temp <- rbind(Plot_lyr_temp[49:96,1:73],Plot_lyr_temp[1:48,1:73])
  Plot_lyr_prec <- rbind(Plot_lyr_prec[49:96,1:73],Plot_lyr_prec[1:48,1:73])
  
  ##### Point Layer
  
  Point_Lyr_temp <- as.tibble(CORR_CAVE) %>% filter(!is.na(p_TEMP)) %>% filter(p_TEMP<0.1) %>% select(entity_id, lon,lat,CORR_TEMP) %>% rename(value = CORR_TEMP)
  Point_Lyr_temp_not <- as.tibble(CORR_CAVE) %>% filter(!entity_id %in% Point_Lyr_temp$entity_id) %>% select(entity_id, lon,lat,CORR_TEMP) %>% rename(value = CORR_TEMP)
  Point_Lyr_prec <- as.tibble(CORR_CAVE) %>% filter(!is.na(p_PREC)) %>% filter(p_PREC<0.1) %>% select(entity_id, lon,lat,CORR_PREC) %>% rename(value = CORR_PREC)
  Point_Lyr_prec_not <- as.tibble(CORR_CAVE) %>% filter(!entity_id %in% Point_Lyr_prec$entity_id) %>% select(entity_id, lon,lat,CORR_PREC) %>% rename(value = CORR_PREC)
  
  Point_Lyr_temp_not_p <- projection_ptlyr(as.data.frame(Point_Lyr_temp_not %>% select(lon,lat,value)), projection = as.character('+proj=robin +datum=WGS84'))
  Point_Lyr_prec_not_p <- projection_ptlyr(as.data.frame(Point_Lyr_prec_not %>% select(lon,lat,value)), projection = as.character('+proj=robin +datum=WGS84'))
  Point_Lyr_temp_p <- projection_ptlyr(as.data.frame(Point_Lyr_temp %>% select(lon,lat,value)), projection = as.character('+proj=robin +datum=WGS84'))
  Point_Lyr_prec_p <- projection_ptlyr(as.data.frame(Point_Lyr_prec %>% select(lon,lat,value)), projection = as.character('+proj=robin +datum=WGS84'))
  
  #Plot
  
  GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE <- 4
  
  NA_plot_lyr = Plot_lyr_temp
  NA_plot_lyr[!is.na(NA_plot_lyr)] = 0
  NA_plot_lyr[is.na(NA_plot_lyr)] = 1
  
  
  plot_temp <- STACYmap_NA(gridlyr = Plot_lyr_temp, centercolor = 0, graticules = T,
                           NA_gridlyr = NA_plot_lyr, NA_color = "grey", legend_names = list(grid = TeX("$\\rho (T, \\delta^{18}O)$"))) +
    geom_point(data = as.data.frame(Point_Lyr_temp_not_p), aes(x = long, y = lat), shape = 20, color = "black", size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE-1.5) +
    geom_point(data = as.data.frame(Point_Lyr_temp_p), aes(x = long, y = lat, fill = layer), shape = 21, color = "black", size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE) +
    theme(panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  plot_temp
  
  NA_plot_lyr = Plot_lyr_prec
  NA_plot_lyr[!is.na(NA_plot_lyr)] = 0
  NA_plot_lyr[is.na(NA_plot_lyr)] = 1
  
  plot_prec <- STACYmap_NA(gridlyr = Plot_lyr_prec, centercolor = 0, graticules = T,
                           NA_gridlyr = NA_plot_lyr, NA_color = "grey",
                           legend_names = list(grid = TeX("$\\rho (P, \\delta^{18}O)$"))) + 
    geom_point(data = as.data.frame(Point_Lyr_prec_not_p), aes(x = long, y = lat), shape = 20, color = "black", size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE-1.5) +
    geom_point(data = as.data.frame(Point_Lyr_prec_p), aes(x = long, y = lat, fill = layer), shape = 21, color = "black", size = GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE) +
    theme(panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 12),
          legend.title = element_text(size = 12))
  
  #plot_prec
  
  library(ggpubr)
  plot <- ggarrange(plot_temp, plot_prec,
                    labels = c("(a)", "(b)"),
                    ncol = 2, nrow = 1)
  
  plot  %>% ggsave(filename = paste0('Fig7_Correlation_xnap',run, '.pdf'), plot = ., path = 'Paper_Plots', 
                   width = 2*12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
  plot  %>% ggsave(filename = paste0('Fig7_Correlation_xnap',run, '.png'), plot = ., path = 'Paper_Plots', 
                   width = 2*12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "png")
}

remove(COR, double_time, plot, Plot_lyr_prec, Plot_lyr_prec_p, Plot_lyr_temp, Plot_lyr_temp_p, plot_prec, plot_temp, Point_Lyr_prec, Point_Lyr_temp)
remove(entity, ii, record, run, sim, var, Point_Lyr_prec_not, Point_Lyr_prec_not_p, Point_Lyr_prec_p, Point_Lyr_temp_not, Point_Lyr_temp_not_p, Point_Lyr_temp_p)
remove(CORR, CORR_CAVE, CORR_FIELD, data_rec, NA_plot_lyr,lon,lat)
