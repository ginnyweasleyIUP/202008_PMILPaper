library(plyr)
library(dplyr)
library(rgdal)

## PLOTTING #####################################

source("Functions/STACYmap_PMIL.R")
source("Functions/STACYmap_PMIL_logscale.R")
source("Functions/aw_mean.R")
source("Functions/projection_ptlyr.R")

GLOBAL_STACY_OPTIONS$GLOBAL_FONT_SIZE = 10

for(run in c("a")){
  ####### DWEQ Plots to double check
  
  CAVElyr <- data.frame(
    lon = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
    lat = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)), 
    value = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
  )
  
  
  for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
    entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
    CAVElyr$lon[ii]   = DATA_past1000$CAVES$entity_info$longitude[ii]
    CAVElyr$lat[ii]   = DATA_past1000$CAVES$entity_info$latitude[ii]
    data = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
    CAVElyr$value[ii] = mean(data[[paste0("d18O_dw_eq_",run)]], na.rm = T)
  }
  
  CAVElyr_used <- data.frame(
    lon = CAVElyr$lon[mask_mean],
    lat = CAVElyr$lat[mask_mean],
    value = CAVElyr$value[mask_mean]
  )
  
  ptlyr <- projection_ptlyr(CAVElyr_used, as.character('+proj=robin +datum=WGS84'))
  
  leg_name = expression(paste(delta^{plain(18)}, plain(O)[plain(rec)], " (%)"))
  
  plot_caves <- STACYmap(coastline = TRUE) +
    geom_point(data = ptlyr, aes(x = long, y = lat, fill = layer), shape = 21, alpha = 0.7, color = "black",
               size = 4, show.legend = c(color =TRUE)) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, 'YlGnBu')), 
                         limits = c(min(ptlyr$layer, na.rm = T), max(ptlyr$layer, na.rm = T)),
                         name = leg_name, 
                         guide = guide_colorbar(barwidth = 10, barheight = 0.3)) +
    theme(legend.direction = "horizontal", 
          panel.border = element_blank(),
          legend.background = element_blank(),
          axis.text = element_blank(),
          text = element_text(size = 10),
          legend.title = element_text(size = 10),
          panel.ontop = F)
  
  plot_caves  %>% ggsave(filename = paste0('Fig3_Mean_CAVESONLY_xnap',run, '.pdf'), plot = ., path = 'CheckUp_Plots', 
                         width = 12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
  plot_caves  %>% ggsave(filename = paste0('Fig3_Mean_CAVESONLY_xnap',run, '.png'), plot = ., path = 'CheckUp_Plots', 
                         width = 12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "png")  
}

rm(plot_caves, ptlyr, leg_name, CAVElyr, CAVElyr_used, entity, data)
