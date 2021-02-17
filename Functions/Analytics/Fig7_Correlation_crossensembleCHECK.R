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

CORR_FIELD <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)),
                   TEMP_P = array(dim = c(96,73)), PREC_P = array(dim = c(96,73)))
# for(var in c("TEMP", "PREC")){
#   for(lon in 1:96){
#     print(lon)
#     for(lat in 2:72){
#       if(sum(is.na(DATA_past1000$SIM_yearly_a$ISOT[lon,lat,]))<100){
#         CORR = cor.test(apply(rbind(DATA_past1000$SIM_yearly_a[[var]][lon,lat,],
#                                     DATA_past1000$SIM_yearly_b[[var]][lon,lat,],
#                                     DATA_past1000$SIM_yearly_c[[var]][lon,lat,]),2,mean, na.rm = T),
#                         apply(rbind(DATA_past1000$SIM_yearly_a$ITPC[lon,lat,],
#                                     DATA_past1000$SIM_yearly_b$ITPC[lon,lat,],
#                                     DATA_past1000$SIM_yearly_c$ITPC[lon,lat,]),2,mean, na.rm = T), na.rm = T)
#         CORR_FIELD[[var]][lon,lat] = CORR$estimate[[1]]
#         CORR_FIELD[[paste0(var, "_P")]][lon,lat] = CORR$p.value
#       }else{
#         CORR_FIELD[[var]][lon,lat] = NA
#         CORR_FIELD[[paste0(var, "_P")]][lon,lat] = NA
#       }
#     }
#   }
# }

for(var in c("TEMP", "PREC")){
  for(lon in 1:96){
    print(lon)
    for(lat in 2:72){
      if(sum(is.na(DATA_past1000$SIM_yearly_a$ISOT[lon,lat,]))<100){
        CORR = cor.test(DATA_past1000$SIM_yearly_a[[var]][lon,lat,], DATA_past1000$SIM_yearly_c$ITPC[lon,lat,], na.rm = T)
        CORR_FIELD[[var]][lon,lat] = CORR$estimate[[1]]
        CORR_FIELD[[paste0(var, "_P")]][lon,lat] = CORR$p.value
      }else{
        CORR_FIELD[[var]][lon,lat] = NA
        CORR_FIELD[[paste0(var, "_P")]][lon,lat] = NA
      }
    }
  }
}


Plot_lyr_temp = CORR_FIELD$TEMP
Plot_lyr_temp_p =CORR_FIELD$TEMP_P
Plot_lyr_temp[Plot_lyr_temp_p > 0.1] = NA
Plot_lyr_prec = CORR_FIELD$PREC
Plot_lyr_prec_p = CORR_FIELD$PREC_P
Plot_lyr_prec[Plot_lyr_prec_p > 0.1] = NA

Plot_lyr_temp <- rbind(Plot_lyr_temp[49:96,1:73],Plot_lyr_temp[1:48,1:73])
Plot_lyr_prec <- rbind(Plot_lyr_prec[49:96,1:73],Plot_lyr_prec[1:48,1:73])

#Plot

GLOBAL_STACY_OPTIONS$GLOBAL_POINT_SIZE <- 4

NA_plot_lyr = Plot_lyr_temp
NA_plot_lyr[!is.na(NA_plot_lyr)] = 0
NA_plot_lyr[is.na(NA_plot_lyr)] = 1


plot_temp <- STACYmap_NA(gridlyr = Plot_lyr_temp, centercolor = 0, graticules = T,
                         NA_gridlyr = NA_plot_lyr, NA_color = "grey", legend_names = list(grid = TeX("$\\rho (T, \\delta^{18}O)$"))) +
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
  theme(panel.border = element_blank(),
        legend.background = element_blank(),
        axis.text = element_blank(),
        text = element_text(size = 12),
        legend.title = element_text(size = 12))

plot_prec

library(ggpubr)
plot <- ggarrange(plot_temp, plot_prec,
                  labels = c("(a)", "(b)"),
                  ncol = 2, nrow = 1)

plot  %>% ggsave(filename = paste0('SF_CHeck_Fig7_Correlation_ac.pdf'), plot = ., path = 'Sup_Plots', 
                 width = 2*12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")
plot  %>% ggsave(filename = paste0('SF_Check_Fig7_Correlation_ac.png'), plot = ., path = 'Sup_Plots', 
                 width = 2*12, height = 12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "png")


remove(COR, double_time, plot, Plot_lyr_prec, Plot_lyr_prec_p, Plot_lyr_temp, Plot_lyr_temp_p, plot_prec, plot_temp, Point_Lyr_prec, Point_Lyr_temp)
remove(entity, ii, record, run, sim, var, Point_Lyr_prec_not, Point_Lyr_prec_not_p, Point_Lyr_prec_p, Point_Lyr_temp_not, Point_Lyr_temp_not_p, Point_Lyr_temp_p)
remove(CORR, CORR_CAVE, CORR_FIELD, data_rec, NA_plot_lyr,lon,lat)
