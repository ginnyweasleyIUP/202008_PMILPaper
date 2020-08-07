#################################################
## Paper Figure SPECTRUM ########################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

#################################################

## PLOT

COLZ <- c("black", "#1A2254", "#0A4296", "#278BCE","#91002B", "#BD6B73", "black")
cex_text = 1

var = "ITPC"

cairo_pdf(file = paste0("Sup_Plots/TS_1-54_",var,".pdf"), width = 21, height = 29.7)
par(mfrow=c(9,6),oma = c(4,3,0,0) + 0.1,mar = c(0,0,1,1) + 0.1)
for(ii in 1:54){
  entity = DATA_past1000$CAVES$entity_info$entity_id[mask_mean][ii]
  data_rec <- DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  yrange = range(c(range(data_rec$d18O_dw_eq_a, na.rm = T), 
                     range(data_rec[[paste0(var, "_a")]], na.rm = T)))
  yrange = c(-15,0)
  if(ii%%6 == 1){
    plot(data_rec$interp_age, data_rec$d18O_measurement, 
         col = adjustcolor(COLZ[1], alpha = 0.5), ylim = yrange, xlim = c(0,1150), type = "l",
         ylab = "",
         xaxt = "n",
         xlab = "", 
         lwd = 2,panel.first = grid())
    mtext("d18O in [‰]", side = 2, line= 2, cex = 0.8)
    
  }else{
    plot(data_rec$interp_age, data_rec$d18O_measurement, 
         col = adjustcolor(COLZ[1], alpha = 0.5), ylim = yrange, xlim = c(0,1150), type = "l",
         ylab = "",
         yaxt = "n",
         xaxt = "n",
         xlab = "", 
         lwd = 2, panel.first = grid())
  }
  
  lines(gaussdetr(zoo(x=data_rec$d18O_measurement, 
                      order.by = data_rec$interp_age), tsc.in = 100)$Xsmooth, col = COLZ[1], lw = 3)
  
  lines(data_rec$interp_age, data_rec$d18O_dw_eq_a, col = adjustcolor(COLZ[2], alpha.f = 0.5), lw = 2)
  lines(gaussdetr(zoo(x = data_rec$d18O_dw_eq_a, 
                      order.by = data_rec$interp_age), tsc.in = 100)$Xsmooth, 
        col = COLZ[3], lw = 2)
  
  lines(data_rec$interp_age, data_rec[[paste0(var, "_a")]], col = adjustcolor(COLZ[5], alpha.f = 0.5), lw = 2)
  lines(gaussdetr(zoo(x = data_rec[[paste0(var, "_a")]], 
                      order.by = data_rec$interp_age), tsc.in = 100)$Xsmooth, 
        col = COLZ[6], lw = 2)
  
  if(ii %in% c(49,50,51,52,53,54)){
    axis(1,at=seq(0,1100,by=200),labels=FALSE,col="black")
    mtext(side=1,at=seq(0,1100,by=200),seq(0,1100,by=200),line = 1.5, cex=cex_text,col="black")
    mtext("years BP", side = 1, line= 3, cex = 0.8)
  }
  
  
  legend("bottomleft", legend = c(paste0("Entity ", entity," d18O"), paste0("Entity ", entity," dweq"), "HadCM3 down-sampled"), 
         col = c(COLZ[1],COLZ[3],COLZ[6]), lwd = c(2,2), lty = c(1,1), bty = "n", cex = cex_text)
  
}



dev.off()

cairo_pdf(file = paste0("Sup_Plots/TS_55-106_",var,".pdf"), width = 21, height = 29.7)
par(mfrow=c(9,6),oma = c(4,3,0,0) + 0.1,mar = c(0,0,1,1) + 0.1)
for(ii in 55:106){
  entity = DATA_past1000$CAVES$entity_info$entity_id[mask_mean][ii]
  data_rec <- DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
  site = DATA_past1000$CAVES$entity_info$site_id[DATA_past1000$CAVES$entity_info$entity_id == entity]
  yrange = range(c(range(data_rec$d18O_dw_eq_a, na.rm = T), 
                   range(data_rec[[paste0(var, "_a")]], na.rm = T)))
  yrange = c(-15,0)
  if(ii%%6 == 1){
    plot(data_rec$interp_age, data_rec$d18O_measurement, 
         col = adjustcolor(COLZ[1], alpha = 0.5), ylim = yrange, xlim = c(0,1150), type = "l",
         ylab = "",
         xaxt = "n",
         xlab = "", 
         lwd = 2,panel.first = grid())
    mtext("d18O in [‰]", side = 2, line= 2, cex = 0.8)
    
  }else{
    plot(data_rec$interp_age, data_rec$d18O_measurement, 
         col = adjustcolor(COLZ[1], alpha = 0.5), ylim = yrange, xlim = c(0,1150), type = "l",
         ylab = "",
         yaxt = "n",
         xaxt = "n",
         xlab = "", 
         lwd = 2, panel.first = grid())
  }
  
  lines(gaussdetr(zoo(x=data_rec$d18O_measurement, 
                      order.by = data_rec$interp_age), tsc.in = 100)$Xsmooth, col = COLZ[1], lw = 3)
  
  lines(data_rec$interp_age, data_rec$d18O_dw_eq_a, col = adjustcolor(COLZ[2], alpha.f = 0.5), lw = 2)
  lines(gaussdetr(zoo(x = data_rec$d18O_dw_eq_a, 
                      order.by = data_rec$interp_age), tsc.in = 100)$Xsmooth, 
          col = COLZ[3], lw = 2)
  
  lines(data_rec$interp_age, data_rec[[paste0(var, "_a")]], col = adjustcolor(COLZ[5], alpha.f = 0.5), lw = 2)
  lines(gaussdetr(zoo(x = data_rec[[paste0(var, "_a")]], 
                      order.by = data_rec$interp_age), tsc.in = 100)$Xsmooth, 
        col = COLZ[6], lw = 2)
  
  if(ii %in% c(101,102,103,104,105,106)){
    axis(1,at=seq(0,1100,by=200),labels=FALSE,col="black")
    mtext(side=1,at=seq(0,1100,by=200),seq(0,1100,by=200),line = 1.5, cex=cex_text,col="black")
    mtext("years BP", side = 1, line= 3, cex = 0.8)
  }
  
  
  legend("bottomleft", legend = c(paste0("Entity ", entity," d18O"), paste0("Entity ", entity," dweq"), "HadCM3 down-sampled"), 
         col = c(COLZ[1],COLZ[3],COLZ[6]), lwd = c(2,2), lty = c(1,1), bty = "n", cex = cex_text)
  
}

dev.off()

rm(entity, site, yrange, data_rec, COLZ, cex_text, var, ii)
