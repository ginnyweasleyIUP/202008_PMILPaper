library(ncdf4)
library(fields)
source("Functions/STACYmap_PMIL.R")

# nc <- nc_open("/home/ginnyweasley/Dokumente/01_Promotion/07_R_Code/202005_Paper/Data/HadCM3/xnapa_pi_0_decade_0.nc")
# lon <- ncvar_get(nc,"longitude")
# lat <- ncvar_get(nc,"latitude")
# pseudo <- ncvar_get(nc,"pseudo")
# surface <- ncvar_get(nc,"surface")
# time <- ncvar_get(nc,"t")
# pftfrac <- ncvar_get(nc,"field1391")
# 
# #HadCM3 Triffid PFTs:
# # 1) Broadleaf tree [Tropical], 2) Needleleaf tree, 3) C3 grass, 4) C4 grass, 5) Shrub, 6) Urban, 7) Inland water, 8) Bare soil, 9) Ice
col_scores = c("#668822", "#117733","#558877", "#88BBAA", "#946846", "#B2B2B2", "#7BAFDE", "#FFDD44","white")
# 
# pft_max <- apply(pftfrac[,,1:9,6],c(1,2),function(x) {if (!is.na(x[1])) {return(which.max(x))} else {NA}})
# pft_max[is.na(pft_max)] = 7 #turn into water
# image.plot(lon,rev(lat),pft_max[,73:1],main="Dominant vegetation",col=col_scores);maps::map(add=T,wrap=c(0,360))

load("Data/LM_HadCM3_JuneVegetation.R")

plot <- STACYmap(gridlyr = rbind(DATA_VEGETATION[49:96,1:73],DATA_VEGETATION[1:48,1:73]),
                 graticules = TRUE,colorscheme = col_scores,
                 legend_cb = F, legend_num_breaks = 9, legend_names = list(grid = "PFT")) +
  theme(panel.border = element_blank(), legend.background = element_blank(), axis.text = element_blank(), legend.text = element_text(size = 12)) 

plot %>% ggsave(filename = paste0('SF_Vegetation_June.pdf'), plot = ., path = 'Sup_Plots', 
                width = 2*12, height = 2*12/8.3*PLOTTING_VARIABLES$HEIGHT, units = 'cm', dpi = 'print', device = "pdf")

#HadCM3 Triffid PFTs:
# 1) Broadleaf tree [Tropical], 2) Needleleaf tree, 3) C3 grass, 4) C4 grass, 5) Shrub, 6) Urban, 7) Inland water, 8) Bare soil, 9) Ice
#pft.names <- c("Broadleaf tree","Needleleaf tree","C3 grass","C4 grass","Shrub","Urban","Inland water","Bare soil","Ice")
rm(DATA_VEGETATION, col_scores)
