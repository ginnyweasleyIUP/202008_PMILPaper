#################################################
## SEASONALITY

library(latex2exp)

# CALCULATION

source("Functions/SubsampleTimeseriesBlock_highresNA.R")

run = "a"

SEASON_CORR <-  list(TEMP = list(entity_id = numeric(sum(mask_var)), lon = numeric(sum(mask_var)), lat = numeric(sum(mask_var)), 
                                 corr = numeric(sum(mask_var)), corr_p = numeric(sum(mask_var)), season = numeric(sum(mask_var))),
                     PREC = list(entity_id = numeric(sum(mask_var)), lon = numeric(sum(mask_var)), lat = numeric(sum(mask_var)), 
                                 corr = numeric(sum(mask_var)), corr_p = numeric(sum(mask_var)), season = numeric(sum(mask_var))), 
                     ISOT = list(entity_id = numeric(sum(mask_var)), lon = numeric(sum(mask_var)), lat = numeric(sum(mask_var)), 
                                 corr = numeric(sum(mask_var)), corr_p = numeric(sum(mask_var)), season = numeric(sum(mask_var))))

for(var in c("TEMP", "PREC", "ISOT")){
  for(ii in 1:sum(mask_var)){
    entity = DATA_past1000$CAVES$entity_info$entity_id[mask_var][ii]
    SEASON_CORR[[var]]$entity_id[ii] = DATA_past1000$CAVES$entity_info$entity_id[mask_var][ii]
    SEASON_CORR[[var]]$lon[ii] = DATA_past1000$CAVES$entity_info$longitude[mask_var][ii]
    SEASON_CORR[[var]]$lat[ii] = DATA_past1000$CAVES$entity_info$latitude[mask_var][ii]
    
    data_rec <- DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
    
    for(season in c("WINTER", "SPRING", "SUMMER", "AUTUMN")){
      
      data_season <- DATA_past1000$CAVES$season_res[[run]][[season]] %>% filter(entity_id == entity)
      TS_season = SubsampleTimeseriesBlock_highresNA(ts(data = rev(data_season[[var]]),
                                            start = LastElement(data_season$year_BP),
                                            end = FirstElement(data_season$year_BP)),
                                         data_rec$interp_age)
      CORR = cor.test(data_rec$d18O_measurement, as.numeric(TS_season))
      if(CORR$p.value<0.1 && abs(CORR$estimate[[1]])>abs(SEASON_CORR[[var]]$corr[ii])){
        SEASON_CORR[[var]]$corr[ii] = CORR$estimate[[1]]
        SEASON_CORR[[var]]$corr_p[ii] = CORR$p.value
        SEASON_CORR[[var]]$season[ii] = season
      }
    }
    #yearly:
    CORR = cor.test(data_rec$d18O_measurement, data_rec[[paste0(var,"_",run)]])
    if(CORR$p.value<0.1 && abs(CORR$estimate[[1]])>abs(SEASON_CORR[[var]]$corr[ii])){
      SEASON_CORR[[var]]$corr[ii] = CORR$estimate[[1]]
      SEASON_CORR[[var]]$corr_p[ii] = CORR$p.value
      SEASON_CORR[[var]]$season[ii] = "YEARLY"
    }
  }
}



#################################################

source("Functions/STACYmap_PMIL.R")
source("Functions/projection_ptlyr.R")

run = "a"


PLOT <- list()
for(var in c("TEMP", "PREC","ISOT")){
  Plot_Lyr <- data.frame(lon = SEASON_CORR[[var]]$lon, 
                         lat = SEASON_CORR[[var]]$lat,
                         layer = SEASON_CORR[[var]]$season,
                         value = SEASON_CORR[[var]]$corr)
  
  Plot_Lyr$layer[Plot_Lyr$layer == "0"] = NA
  
  mask_china <- logical(length(Plot_Lyr$lon))
  
  for(ii in 1:length(Plot_Lyr$lon)){
    if(is.na(Plot_Lyr$lon[ii])){next}
    if(Plot_Lyr$lon[ii] > 100 & Plot_Lyr$lon[ii] < 120){
      if(Plot_Lyr$lat[ii] < 35 & Plot_Lyr$lat[ii] > 22){
        mask_china[ii] = T}
    }
  }
  
  ptlyr_china <- data.frame(
    lon = Plot_Lyr$lon[mask_china],
    lat = Plot_Lyr$lat[mask_china],
    layer = Plot_Lyr$layer[mask_china],
    value = Plot_Lyr$value[mask_china]
  )
  
  ptlyr_rest <- data.frame(
    lon = Plot_Lyr$lon[!mask_china],
    lat = Plot_Lyr$lat[!mask_china],
    layer = Plot_Lyr$layer[!mask_china],
    value = Plot_Lyr$value[!mask_china]
  )
  
  
  ptlyr_china_p <- projection_ptlyr(ptlyr_china, as.character('+proj=robin +datum=WGS84'))#
  ptlyr_rest_p <- projection_ptlyr(ptlyr_rest, as.character('+proj=robin +datum=WGS84'))
  ptlyr_china_p$value <- abs(ptlyr_china$value)
  ptlyr_rest_p$value <- abs(ptlyr_rest$value)
  
  ptlyr_china_p$layer <- factor(ptlyr_china$layer)
  ptlyr_rest_p$layer <- factor(ptlyr_rest$layer)
  
  if(var =="TEMP"){title = "Temperature"}
  if(var =="PREC"){title = "Precipitation"}
  if(var =="ISOT"){title = "Isotopic composition"}
  
  
  plot <- STACYmap(coastline = T) +
    geom_point(data = ptlyr_china_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
               show.legend = c(color =TRUE, size = TRUE), position = position_jitter(width = 1000000, height = 500000)) +
    geom_point(data = ptlyr_rest_p, aes(x = long, y = lat, color = layer, size = value), shape = 19, 
               show.legend = c(color =TRUE, size = TRUE)) +
    scale_color_manual(name = "Seasons", labels = c("DJF", "MAM", "JJA", "SON", "year"), 
                       
                       values = c("#0072b2", "#009e73", "#d55e00", "#f0e442", "#483D3F")) +
    scale_size(name = "abs. corr.") +
    ggtitle(title) +
    theme(panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          legend.direction = "horizontal",
          text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5))+
    guides(color=guide_legend(nrow=2,byrow=TRUE), size = guide_legend(nrow = 2, byrow = T))
  
  #plot
  # plot  %>% ggsave(filename = paste('Paper_Plot_A2_Seasons_prec', 'pdf', sep = '.'), plot = ., path = 'Plots', 
  #                  width = 2*8.3, height =2*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
  PLOT[[var]] <- plot
}

library(ggpubr)
plot <- ggarrange(PLOT$TEMP, PLOT$PREC, PLOT$ISOT,
                  labels = c("(a)", "(b)", "(c)"),
                  ncol = 1, nrow = 3)

plot  %>% ggsave(filename = paste0('SF_Map1_Seasons_',run, '.pdf'), plot = ., path = 'Sup_Plots', 
                 width = 2*PLOTTING_VARIABLES$WIDTH, height = 4*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
plot  %>% ggsave(filename = paste0('SF_Map1_Seasons_',run, '.png'), plot = ., path = 'Sup_Plots', 
                 width = 2*PLOTTING_VARIABLES$WIDTH, height = 4*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "png")

#################################################
##CORR Map

load("Data/LM_HadCM3_correlation.R")

CORR_FIELD <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)),
                   TEMP_P = array(dim = c(96,73)), PREC_P = array(dim = c(96,73)))

CORR_FIELD$TEMP = DATA_EXPORT_CORR[[run]]$TEMP
CORR_FIELD$TEMP_p = DATA_EXPORT_CORR[[run]]$TEMP_P
CORR_FIELD$PREC = DATA_EXPORT_CORR[[run]]$PREC
CORR_FIELD$PREC_p = DATA_EXPORT_CORR[[run]]$PREC_P

rm(DATA_EXPORT_CORR)

CORR_CAVE <- list(entity_id = SEASON_CORR$TEMP$entity_id,
                  lon = SEASON_CORR$TEMP$lon,
                  lat = SEASON_CORR$TEMP$lat, 
                  value_temp = SEASON_CORR$TEMP$corr,
                  value_prec = SEASON_CORR$PREC$corr)

CORR_CAVE$value_temp[SEASON_CORR$TEMP$season == "0"] = NA
CORR_CAVE$value_prec[SEASON_CORR$PREC$season == "0"] = NA


#################################################
## PLOTS ########################################
#################################################

Plot_lyr_temp = CORR_FIELD$TEMP
Plot_lyr_temp_p =CORR_FIELD$TEMP_P
Plot_lyr_temp[Plot_lyr_temp_p > 0.1] = NA
Plot_lyr_prec = CORR_FIELD$PREC
Plot_lyr_prec_p = CORR_FIELD$PREC_P
Plot_lyr_prec[Plot_lyr_prec_p > 0.1] = NA
  
Plot_lyr_temp <- rbind(Plot_lyr_temp[49:96,1:73],Plot_lyr_temp[1:48,1:73])
Plot_lyr_prec <- rbind(Plot_lyr_prec[49:96,1:73],Plot_lyr_prec[1:48,1:73])
  
##### Point Layer

Point_Lyr_temp <- as.tibble(CORR_CAVE) %>% select(lon,lat,value_temp) %>% filter(!is.na(value_temp)) %>% rename(value = value_temp)
Point_Lyr_temp_not <- as.tibble(CORR_CAVE)  %>% select(lon,lat, value_temp) %>% filter(is.na(value_temp)) %>% rename(value = value_temp)
Point_Lyr_prec <- as.tibble(CORR_CAVE) %>% select(lon,lat,value_prec) %>% filter(!is.na(value_prec)) %>% rename(value = value_prec)
Point_Lyr_prec_not <- as.tibble(CORR_CAVE)  %>% select(lon,lat, value_prec) %>% filter(is.na(value_prec)) %>% rename(value = value_prec)


Point_Lyr_temp_not_p <- projection_ptlyr(as.data.frame(Point_Lyr_temp_not %>% select(lon,lat,value)), projection = as.character('+proj=robin +datum=WGS84'))
Point_Lyr_prec_not_p <- projection_ptlyr(as.data.frame(Point_Lyr_prec_not %>% select(lon,lat,value)), projection = as.character('+proj=robin +datum=WGS84'))
Point_Lyr_temp_p <- projection_ptlyr(as.data.frame(Point_Lyr_temp %>% select(lon,lat,value)), projection = as.character('+proj=robin +datum=WGS84'))
Point_Lyr_prec_p <- projection_ptlyr(as.data.frame(Point_Lyr_prec %>% select(lon,lat,value)), projection = as.character('+proj=robin +datum=WGS84'))

#Plot

source("Functions/STACYmap_PMIL_NAgrid.R")  
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

#plot_temp

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

plot  %>% ggsave(filename = paste0("SF_Map2_Seasons_",run,".pdf"), plot = ., path = 'Sup_Plots', 
                 width = 2*12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
plot  %>% ggsave(filename = paste0("SF_Map2_Seasons_",run,".png"), plot = ., path = 'Sup_Plots', 
                 width = 2*12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "png")

remove(plot, Plot_lyr_prec, Plot_lyr_prec_p, Plot_lyr_temp, Plot_lyr_temp_p, plot_prec, plot_temp, Point_Lyr_prec, Point_Lyr_temp)
remove(Point_Lyr_prec_not, Point_Lyr_prec_not_p, Point_Lyr_prec_p, Point_Lyr_temp_not, Point_Lyr_temp_not_p, Point_Lyr_temp_p)
remove(CORR_CAVE, CORR_FIELD, NA_plot_lyr, Plot_Lyr, ptlyr_china,  ptlyr_china_p, ptlyr_rest,ptlyr_rest_p)
rm(entity, ii, mask_china, run, season, title, TS_season, var, PLOT, data_rec, data_season, CORR)
