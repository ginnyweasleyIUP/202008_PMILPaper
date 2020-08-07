#################################################
## Paper Figure 8 ###############################
#################################################

## Network 

#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(nest)
library(PaleoSpec)

lats = DATA_past1000$CAVES$entity_info$latitude[mask_spec]
longs = DATA_past1000$CAVES$entity_info$longitude[mask_spec]

dist<-fossil::earth.dist(cbind(lats,longs),dist=TRUE)
dist_matrix <- as.matrix(dist)

#################################################

source("Functions/networkmap_simple3.R")

run = "a"

for(plot in c("full", "down")){
  
  plot_dist = dist_matrix
  plot_dist[lower.tri(dist_matrix)] = NA
  plot_dist_e = dist_matrix_e
  plot_dist_e[lower.tri(dist_matrix_e)] = NA
  
  link_density = 0.05
  C_SIM = NETWORK[[paste0("sim_", plot, "_", run)]]
  C_REC = NETWORK$record
  
  o_sim = order(abs(C_SIM), na.last = F)
  o_rec = order(abs(C_REC), na.last = F)
  C_SIM[o_sim[1:floor((length(o_sim)-link_density*length(o_sim)))]] =  NA
  C_REC[o_rec[1:floor((length(o_rec)-link_density*length(o_rec)))]] = NA
  
  C_SIM[NETWORK[[paste0("sim_", plot, "_", run, "_p")]]>0.1] = NA
  C_REC[NETWORK$record_p>0.1] = NA
  
  plot_c_sim = NETWORK[[paste0("sim_", plot, "_", run)]]
  plot_c_sim[lower.tri(NETWORK[[paste0("sim_", plot, "_", run)]], diag = FALSE)] = NA

  plot_c_rec = NETWORK$record
  plot_c_rec[lower.tri(NETWORK$record)] = NA
  
  plot_c_rec_max_e <- c_ensemble
  plot_c_rec_max_e[lower.tri(c_ensemble)] = NA
  
  boxes_sim <- list()
  boxes_rec <- list()
  
  for(ii in 1:20){
    boxes_sim[[paste0(ii*1000)]] <- na.omit(as.numeric(plot_c_sim[plot_dist<ii*1000 & plot_dist>(ii-1)*1000]))
    boxes_rec[[paste0(ii*1000)]] <- na.omit(as.numeric(plot_c_rec[plot_dist<ii*1000 & plot_dist>(ii-1)*1000]))
  }
  
  all_scale = 1.4
  scaling = 1.5*all_scale
  spacing = 0.7*all_scale
  namcex = 1*all_scale
  colorbar_length = 0.7
  
  # vielleicht erst die stärksten Links aussuchen udn dann unsignifikante raus schmeißen! Damit die gleiche Anzahl an Links da ist!
  rbPal <- colorRampPalette(c("#2166ac", "grey", "#b2182b"))
  COLZ <- array(rbPal(9))
  
  
  pdf(file = paste0("Paper_Plots/Fig8_Network_",plot, "_", run,".pdf"), height= PLOTTING_VARIABLES$HEIGHT, width = PLOTTING_VARIABLES$WIDTH)
  #png(file = paste0("Plots/Paper_Plot_6_Network_a_ds_xnap",run,".png"), height= 100*PLOTTING_VARIABLES$HEIGHT, width = 100*PLOTTING_VARIABLES$WIDTH)
  par(mfrow=c(2,2), mai = c(rep(spacing, 4)), mar = c(1.5,2.5,0,0.5), oma = c(2,0.5,0.5,0.5))
  #SIM MAP
  networkmap_simple3(CMAT = C_SIM, 
                     lat = lats, 
                     lon = longs, 
                     title = "",
                     thresh = 0.1)
  fields::colorbar.plot(x = 190,y = 10, col = rev(COLZ), 
                        strip = c(1,0.75, 0.5, 0.25,0,-0.25, -0.5, -0.75,-1), horizontal = F, strip.length = colorbar_length, strip.width = 0.04)
  axis(4,at=seq(80,-60,by=-2*140/8),labels=FALSE)
  mtext("HadCM3@SISALv2", side = 1, cex = namcex, line = -2, font = 2)
  mtext("(a)", side = 3, adj = 0, cex = namcex, line = -1.5, at = -185)
  #SIM Cor-Dist
  plot(plot_dist[seq(1,length(plot_dist), by = 2)], plot_c_sim[seq(1,length(plot_c_sim), by = 2)], 
       ylim = c(-1,1),xlim = c(0,20000),
       ylab = "", xlab = "", 
       cex = 1, lwd = 0.5, pch = 23,
       panel.first = grid(), col = adjustcolor("grey", alpha.f = 0.3), xaxt = "n")
  limits = c(1, 0.875, 0.625, 0.375, 0.125, -0.125, -0.375, -0.625, -0.875)
  for(ii in 1:(length(limits)-1)){
    points(dist_matrix[C_SIM>limits[ii+1] & C_SIM<limits[ii]], C_SIM[C_SIM>limits[ii+1] & C_SIM<limits[ii]], pch = 23, 
           col = adjustcolor(rev(COLZ)[ii], alpha.f = 0.6))  
  }
  points(dist_matrix[C_SIM<limits[9]], C_SIM[C_SIM<limits[9]], pch = 23, col = adjustcolor(rev(COLZ)[9], alpha.f = 0.7))
  abline(h=0)
  for(ii in 1:20){
    boxplot(boxes_sim[[paste0(ii*1000)]], add = TRUE, at = c(ii*1000-500),boxwex = 1000, names = "n", axes = F, outline = F)  
  }
  
  
  lo <- loess(plot_c_sim[order(plot_dist)] ~ plot_dist[order(plot_dist)], span = 0.2)

  lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
  
  #lines(lowess(lowess_dist_sorted,lowess_c_sim_sorted, f=0.1), lwd = 4, col = "#B2182B")
  mtext("(b)", side = 3, adj = 0, cex = namcex, line = -1.5, at = 1000)
  mtext("HadCM3@SISALv2", side = 1, cex = namcex, line = -2, font = 2)
  
  #SISAL MAP
  networkmap_simple3(CMAT = C_REC, 
                     lat = lats, 
                     lon = longs,
                     title = "",
                     thresh = 0.1)
  fields::colorbar.plot(x = 190,y = 10, col = rev(COLZ), 
                        strip = c(1,0.75, 0.5, 0.25,0,-0.25, -0.5, -0.75,-1), horizontal = F, strip.length = colorbar_length, strip.width = 0.04)
  axis(4,at=seq(80,-60,by=-2*140/8),labels=FALSE)
  mtext("SISALv2", side = 1, cex = namcex, line = -2, font = 2)
  mtext("(c)", side = 3, adj = 0, cex = namcex, line = -1.5, at = -185)
  
  #SISAL Cor-Dist
  plot(plot_dist[seq(1,length(plot_dist), by = 2)], plot_c_rec[seq(1,length(plot_c_rec), by = 2)], 
       ylim = c(-1,1),
       xlim = c(0,20000),
       ylab = "",
       xlab = "",
       cex.axis = all_scale, 
       lwd = 1,  pch = 24,
       panel.first = grid(), col = adjustcolor("grey", alpha.f = 0.3))
  limits = c(1, 0.875, 0.625, 0.375, 0.125, -0.125, -0.375, -0.625, -0.875)
  for(ii in 1:(length(limits)-1)){
    points(dist_matrix[C_REC>limits[ii+1] & C_REC<=limits[ii]], C_REC[C_REC>limits[ii+1] & C_REC<=limits[ii]], pch = 24, 
           col = adjustcolor(rev(COLZ)[ii], alpha.f = 0.6))  
  }
  points(dist_matrix[C_REC<limits[9]], C_REC[C_REC<limits[9]], pch = 24, col = adjustcolor(rev(COLZ)[9], alpha.f = 0.7))
  abline(h=0)
  for(ii in 1:20){
    boxplot(boxes_rec[[paste0(ii*1000)]], add = TRUE, at = c(ii*1000-500),boxwex = 1000, names = "n", axes = F, outline = F)  
  }
  
  lo <- loess(plot_c_rec[order(plot_dist)] ~ plot_dist[order(plot_dist)], span = 0.2)
  lines(lo$x, lo$fitted, lwd = 4, col = "#B2182B")
  lo <- loess(plot_c_rec_max_e[order(plot_dist_e)] ~ plot_dist_e[order(plot_dist_e)], span = 0.2)
  lines(lo$x, lo$fitted, lwd = 4, col = "#0068C4")
  mtext("SISALv2", side = 1, cex = namcex, line = -2, font = 2)
  mtext("Distance between pairs (km)", side= 1, line = 2)
  mtext("(d)", side = 3, adj = 0, cex = namcex, line = -1.5, at = 0)
  text(20000, 0.9, "original chron.", col = "#B2182B", adj = 1, cex = all_scale)
  text(20000, 0.77, "sisal ensemble", col = "#0068C4", adj = 1, cex = all_scale)
  
  dev.off()
  
}




rm(boxes_rec, boxes_sim, COLZ, lo, plot_c_rec, plot_dist, ii, lats, link_density,longs, namcex, o_rec, o_sim, plot, scaling, spacing,plot_c_sim,
   C_REC, C_SIM, plot_c_rec_max_e, plot_dist_e, colorbar_length, limits, all_scale, dist, run)
