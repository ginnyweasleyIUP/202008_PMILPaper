#################################################
## Paper Figure 3 ###############################
#################################################

library(plyr)
library(dplyr)
library(rgdal)

## PLOTTING #####################################

source("Functions/STACYmap_PMIL.R")
source("Functions/STACYmap_PMIL_logscale.R")
source("Functions/aw_mean.R")
source("Functions/projection_ptlyr.R")

GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 12
txt_size = 12

load("Data/LM_HadCM3_annualmean.R")

for(run in c("a")){
  for(var in c("ISOT", "ITPC")){
    # TEMP LAYER
    temp_lyr <- DATA_EXPORT_FIELD[[run]]$TEMP_MEAN
    plot_temp <- STACYmap(gridlyr = rbind(temp_lyr[49:96,1:73],temp_lyr[1:48,1:73]),
                          legend_names = list(grid = "Temperature (°C)"),graticules = TRUE,colorscheme = "temp", centercolor = 0) +
      theme(panel.border = element_blank(), legend.background = element_blank(), axis.text = element_blank(), legend.text = element_text(size = txt_size)) 
    
    # PREC LAYER
    prec_lyr <- DATA_EXPORT_FIELD[[run]]$PREC_MEAN
    plot_prec <- STACYmap(gridlyr = 8.6148e4*rbind(prec_lyr[49:96,1:73],prec_lyr[1:48,1:73]),
                          legend_names = list(grid = "Precipitation (mm/day)"), graticules = TRUE, colorscheme = "prcp_grd") +
      theme(panel.border = element_blank(), legend.background = element_blank(), axis.text = element_blank(), legend.text = element_text(size = txt_size)) 
    
    
    # ISOT POINTS
    CAVElyr_isot <- data.frame(
      lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
      lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)), 
      value = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
    )
    
    
    for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
      entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
      CAVElyr_isot$lon[ii]   = DATA_past1000$CAVES$entity_info$longitude[ii]
      CAVElyr_isot$lat[ii]   = DATA_past1000$CAVES$entity_info$latitude[ii]
      data = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
      CAVElyr_isot$value[ii] = mean(data[[paste0("d18O_dw_eq_",run)]], na.rm = T)
    }
    
    CAVElyr_isot_used <- data.frame(
      lon = CAVElyr_isot$lon[mask_mean],
      lat = CAVElyr_isot$lat[mask_mean],
      value = CAVElyr_isot$value[mask_mean]
    )
    
    remove(entity, ii)
    
    #ISOT LAYER
    isot_lyr <- DATA_EXPORT_FIELD[[run]][[paste0(var,"_MEAN")]]
    Plot_lyr1 <-rbind(isot_lyr[49:96,1:73],isot_lyr[1:48,1:73])
    
    Plot_lyr3 <- Plot_lyr1
    Plot_lyr3[is.na(Plot_lyr3)] = 1000
    Plot_lyr3[Plot_lyr3>0] <- Plot_lyr3[Plot_lyr3>0]+1 
    Plot_lyr3[Plot_lyr3<0] <- Plot_lyr3[Plot_lyr3<0]-1
    Plot_lyr3[Plot_lyr3>0] <- log10(Plot_lyr3[Plot_lyr3>0])
    Plot_lyr3[Plot_lyr3<0] <- - log10(abs(Plot_lyr3[Plot_lyr3<0]))
    Plot_lyr3[abs(Plot_lyr3)>5] <- NA
    Plot_lyr3[,1] <- NA
    Plot_lyr3[,73] <- NA
    
    
    Point_Lyr <- data.frame(
      lon = CAVElyr_isot_used$lon,
      lat = CAVElyr_isot_used$lat,
      value = - log10(abs(CAVElyr_isot_used$value -1))
    )
    
    #Point_Lyr$value[[57]] <- log10(CAVElyr_isot_used$value[[57]]+1)
    
    
    allmax = - min(Plot_lyr3, na.rm = T)
    allmax_real = - min(Plot_lyr1, na.rm = T)
    
    GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE = 4
    
    leg_name = expression(paste(delta^{plain(18)}, plain(O), " (%)"))
    if(var == "ITPC"){leg_name = expression(paste(delta^{plain(18)}, plain(O)[plain(pw)], " (%)"))}
    
    plot_isot <- STACYmap_isot(gridlyr = Plot_lyr3,
                               ptlyr = Point_Lyr,
                               legend_names = list(grid =leg_name),
                               graticules = TRUE,
                               colorscheme = rev(RColorBrewer::brewer.pal(11, 'BrBG'))[c(1:8,10)],
                               #colorscheme = PLOTTING_VARIABLES$COLORS$SLPR,
                               breaks_isot = c(-log10(71),-log10(11),-log10(6), -log10(2), 0, log10(2), log10(6)),
                               labels_isot = c(-70, -10,-5, -1, 0, 1, 5),
                               min_grid = -log10(71), max_grid = log10(12)) +
      theme(panel.border = element_blank(),
            legend.background = element_blank(),
            axis.text = element_blank(),
            text = element_text(size = txt_size)) 
    
    plot_isot
    
    # PLOT MEAN BIAS
    
    CAVElyr_diff <- data.frame(
      lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
      lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)), 
      value = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
    )
    
    
    for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
      entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
      CAVElyr_diff$lon[ii]   = DATA_past1000$CAVES$entity_info$longitude[ii]
      CAVElyr_diff$lat[ii]   = DATA_past1000$CAVES$entity_info$latitude[ii]
      data = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
      CAVElyr_diff$value[ii] = mean(data[[paste0(var,"_",run)]] - data[[paste0("d18O_dw_eq_",run)]], na.rm = T)
    }
    
    CAVElyr_diff_used <- data.frame(
      lon = CAVElyr_diff$lon[mask_mean],
      lat = CAVElyr_diff$lat[mask_mean],
      value = CAVElyr_diff$value[mask_mean]
    )
    #PLOTTING_DATA$FIG3$CAVElyr_diff_used$value[57] = NA
    
    ptlyr <- projection_ptlyr(CAVElyr_diff_used, as.character('+proj=robin +datum=WGS84'))
    
    leg_name = expression(paste(delta^{plain(18)}, plain(O)[plain(sim)]," - ", delta^{plain(18)}, plain(O)[plain(rec)], " (%)"))
    if(var == "ITPC"){leg_name = expression(paste(delta^{plain(18)}, plain(O)["sim,pw"]," - ", delta^{plain(18)}, plain(O)[plain(rec)], " (%)"))}
    
    plot_diff <- STACYmap(coastline = TRUE) +
      geom_point(data = ptlyr, aes(x = long, y = lat, fill = layer), shape = 21, alpha = 0.8, color = "black",
                 size = 4, show.legend = c(color =TRUE)) +
      scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, 'RdBu')), 
                           limits = c(-max(abs(ptlyr$layer), na.rm = T), max(abs(ptlyr$layer), na.rm = T)),
                           name = leg_name, 
                           guide = guide_colorbar(barwidth = 10, barheight = 0.3)) +
      theme(legend.direction = "horizontal", 
            panel.border = element_blank(),
            legend.background = element_blank(),
            axis.text = element_blank(),
            text = element_text(size = txt_size),
            legend.title = element_text(size = txt_size),
            panel.ontop = F)
    
    plot_diff

    
    library(ggpubr)
    plot <- ggarrange(plot_temp, plot_prec, plot_isot, plot_diff,
                      labels = c("(a)", "(b)", "(c)", "(d)"),
                      ncol = 2, nrow = 2)
    
    plot  %>% ggsave(filename = paste0('Fig3_Mean_',var,'_xnap',run, '.pdf'), plot = ., path = 'Paper_Plots', 
                     width = 2*12, height = 2*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
    plot  %>% ggsave(filename = paste0('Fig3_Mean_',var,'_xnap',run, '.png'), plot = ., path = 'Paper_Plots', 
                     width = 2*12, height = 2*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "png")
    
  
    
  }
}

# temp <- simpleawmean(DATA_past1000$SIM_mean$TEMP, seq(from = -90, to = 90, length.out = 73))
# temp_sd <- simpleawsd(DATA_past1000$SIM_mean$TEMP, seq(from = -90, to = 90, length.out = 73))
# prec <- simpleawmean(DATA_past1000$SIM_mean$PREC, seq(from = -90, to = 90, length.out = 73))*8.6148e4
# prec_sd <- simpleawsd(DATA_past1000$SIM_mean$PREC, seq(from = -90, to = 90, length.out = 73))*8.6148e4
# 
# isot <- simpleawmean(DATA_past1000$SIM_mean$ISOT, seq(from = -90, to = 90, length.out = 73))
# isot_sd <- simpleawsd(DATA_past1000$SIM_mean$ISOT, seq(from = -90, to = 90, length.out = 73))
# slpr <- simpleawmean(DATA_past1000$SIM_mean$SLPR, seq(from = -90, to = 90, length.out = 73))
# slpr_sd <- simpleawsd(DATA_past1000$SIM_mean$SLPR, seq(from = -90, to = 90, length.out = 73))




remove(plot_temp, plot_prec, plot_isot, plot_diff, leg_name, Plot_lyr1, Plot_lyr3, plot, Point_Lyr)
remove(temp_lyr, prec_lyr, ptlyr, CAVElyr_diff, CAVElyr_diff_used, CAVElyr_isot, CAVElyr_isot_used)
remove(data, isot_lyr, allmax, allmax_real, entity, ii, run, var, DATA_EXPORT_FIELD)
