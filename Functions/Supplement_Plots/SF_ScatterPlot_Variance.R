#################################################
## MEAN SCATTER PLOTS ###########################
#################################################


library(latex2exp)

No.digits = 2

scatter_data = array(dim = c(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var]),16))
colnames(scatter_data) = c("entity_id", "site_id", "vr_full", "vr_down",
                           "elevation", "latitude", "mean_temp", "mean_prec", "winter_mean_prec", "summer_mean_prec", "elevation_diff", 
                           "dist_entrance", "geology", "cover_thickness", "d18Oc", "mineralogy")

scatter_data[,1] = DATA_past1000$CAVES$entity_info$entity_id[mask_var]
for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])){
  entity = scatter_data[[ii,1]]
  scatter_data[ii,2] = DATA_past1000$CAVES$entity_info$site_id[mask_var][ii]
  data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
  data_yearly_a = DATA_past1000$CAVES$yearly_res$a %>% filter(entity_id == entity)
  data_yearly_b = DATA_past1000$CAVES$yearly_res$b %>% filter(entity_id == entity)
  data_yearly_c = DATA_past1000$CAVES$yearly_res$c %>% filter(entity_id == entity)
  data_seasonal_a = DATA_past1000$CAVES$season_res$a$SUMMER %>% filter(entity_id == entity)
  data_seasonal_b = DATA_past1000$CAVES$season_res$b$SUMMER %>% filter(entity_id == entity)
  data_seasonal_c = DATA_past1000$CAVES$season_res$c$SUMMER %>% filter(entity_id == entity)
  
  scatter_data[ii,3] = as.numeric(mean(c(var(data_rec$d18O_dw_eq_a, na.rm = T)/var(data_yearly_a$ITPC, na.rm = T),
                                         var(data_rec$d18O_dw_eq_b, na.rm = T)/var(data_yearly_b$ITPC, na.rm = T),
                                         var(data_rec$d18O_dw_eq_c, na.rm = T)/var(data_yearly_c$ITPC, na.rm = T)), na.rm = T))
  scatter_data[ii,4] = as.numeric(mean(c(var(data_rec$d18O_dw_eq_a, na.rm = T)/var(data_rec$ITPC_a, na.rm = T),
                                         var(data_rec$d18O_dw_eq_b, na.rm = T)/var(data_rec$ITPC_b, na.rm = T),
                                         var(data_rec$d18O_dw_eq_c, na.rm = T)/var(data_rec$ITPC_c, na.rm = T)), na.rm = T))
  scatter_data[ii,5] = as.numeric(DATA_past1000$CAVES$entity_info$elevation[DATA_past1000$CAVES$entity_info$entity_id == entity])
  scatter_data[ii,6] = as.numeric(DATA_past1000$CAVES$entity_info$latitude[DATA_past1000$CAVES$entity_info$entity_id == entity])
  scatter_data[ii,7] = as.numeric(mean(c(data_rec$TEMP_a, data_rec$TEMP_b, data_rec$TEMP_c), na.rm = T))
  scatter_data[ii,8] = as.numeric(mean(c(data_rec$PREC_a, data_rec$PREC_b, data_rec$PREC_c), na.rm = T))
  scatter_data[ii,9] = as.numeric(DATA_past1000$CAVES$entity_info$winter_prec[DATA_past1000$CAVES$entity_info$entity_id == entity])
  scatter_data[ii,10] = as.numeric(mean(c(data_seasonal_a$PREC, data_seasonal_b$PREC, data_seasonal_c$PREC), na.rm = T))
  scatter_data[ii,11] = as.numeric(DATA_past1000$CAVES$entity_info$elevation_sim[DATA_past1000$CAVES$entity_info$entity_id == entity]-
                                    DATA_past1000$CAVES$entity_info$elevation[DATA_past1000$CAVES$entity_info$entity_id == entity])
  scatter_data[ii,12] = as.numeric(DATA_past1000$CAVES$entity_info$distance_entrance[DATA_past1000$CAVES$entity_info$entity_id == entity])
  scatter_data[ii,13] = as.character(DATA_past1000$CAVES$entity_info$geology[DATA_past1000$CAVES$entity_info$entity_id == entity])
  scatter_data[ii,14] = as.numeric(DATA_past1000$CAVES$entity_info$cover_thickness[DATA_past1000$CAVES$entity_info$entity_id == entity])
  scatter_data[ii,15] = as.numeric(var(data_rec$d18O_measurement, na.rm = T))
  scatter_data[ii,16] = as.character(DATA_past1000$CAVES$entity_info$mineralogy[DATA_past1000$CAVES$entity_info$entity_id == entity])
}

scatter_data = as.data.frame(scatter_data)
scatter_data$entity_id = as.numeric(levels(scatter_data$entity_id))[scatter_data$entity_id]
scatter_data$site_id = as.numeric(levels(scatter_data$site_id))[scatter_data$site_id]
scatter_data$vr_full = as.numeric(levels(scatter_data$vr_full))[scatter_data$vr_full]
scatter_data$vr_down = as.numeric(levels(scatter_data$vr_down))[scatter_data$vr_down]
scatter_data$elevation = as.numeric(levels(scatter_data$elevation))[scatter_data$elevation]
scatter_data$latitude = as.numeric(levels(scatter_data$latitude))[scatter_data$latitude]
scatter_data$mean_temp = as.numeric(levels(scatter_data$mean_temp))[scatter_data$mean_temp]
scatter_data$mean_prec = as.numeric(levels(scatter_data$mean_prec))[scatter_data$mean_prec]*8.6148e4
scatter_data$winter_mean_prec = as.numeric(levels(scatter_data$winter_mean_prec))[scatter_data$winter_mean_prec]
scatter_data$summer_mean_prec = as.numeric(levels(scatter_data$summer_mean_prec))[scatter_data$summer_mean_prec]*8.6148e4
scatter_data$elevation_diff = as.numeric(levels(scatter_data$elevation_diff))[scatter_data$elevation_diff]
scatter_data$dist_entrance = as.numeric(levels(scatter_data$dist_entrance))[scatter_data$dist_entrance]
scatter_data$geology = as.character(levels(scatter_data$geology))[scatter_data$geology]
scatter_data$cover_thickness = as.numeric(levels(scatter_data$cover_thickness))[scatter_data$cover_thickness]
scatter_data$d18Oc = as.numeric(levels(scatter_data$d18Oc))[scatter_data$d18Oc]
scatter_data$mineralogy = as.character(levels(scatter_data$mineralogy))[scatter_data$mineralogy]

# mask for aragonite and calcite
mask_mean_calcite = logical(length = length(scatter_data$entity_id))
mask_mean_aragonite  = logical(length = length(scatter_data$entity_id))

for(ii in 1:length(scatter_data$entity_id)){
  entity = scatter_data$entity_id[ii]
  if(DATA_past1000$CAVES$entity_info$mineralogy[DATA_past1000$CAVES$entity_info$entity_id == entity] == "calcite") {mask_mean_calcite[ii] = T}
  else{mask_mean_aragonite[ii] = T}
}

cexaxis = 1.5

##diff_down
cairo_pdf(width=8,height=10,file= "Sup_Plots/SF_ScatterMean_vr_down.pdf")
#png(file = "Plots/Appendix/A1_ScatterMean_vr_down.png", width = 50*8, height = 50*10)
par(mfrow=c(4,3),oma = c(1,4,0,0) + 0.1,mar = c(3,0,1,1) + 0.1)
#latitude d18Oc
plot(scatter_data$latitude[mask_mean_calcite], scatter_data$d18Oc[mask_mean_calcite],
     xlab = "", ylab = "", ylim = c(0,2), panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
points(scatter_data$latitude[mask_mean_aragonite], scatter_data$d18Oc[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
lm.out <- lm(scatter_data$d18Oc ~ scatter_data$latitude)
abline(lm.out, lwd = 4, col = "#B2182B")
#lines(lowess(scatter_data$latitude[order(scatter_data$latitude)], scatter_data$d18Oc[order(scatter_data$latitude)], f = 2/3, delta = 0.01*diff(range(scatter_data$latitude, na.rm = T))), lwd = 4, 
#      col = adjustcolor("#B2182B", alpha.f = 0.5))
mtext(text = "latitude (°N)",side = 1,line = 2)
mtext(text = TeX("$var(\\delta^{18}O_{calcite})$ (‰$^2$)"),side = 2,line = 2)
text(10, 1, "calcite", cex = 1.5)
text(-10, 1., "aragonite", col = "blue", cex = 1.5)
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
mtext(text = "(a)", side = 3, line = -1.5, adj = 0.02, cex = 1)


#temperature d18Oc
plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$d18Oc[mask_mean_calcite],
     yaxt = 'n', xlab = "", ylab = "", ylim = c(0,2), panel.first = grid(),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$d18Oc[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
lm.out <- lm(scatter_data$d18Oc ~ scatter_data$mean_temp)
abline(lm.out, lwd = 4, col = "#B2182B")
mtext(text = "Mean temperature (°C)",side = 1,line = 2)
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
mtext(text = "(b)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#empty plot
plot(1,1, axes=FALSE,type="n",xlab="",ylab="")

#latitude
plot(scatter_data$latitude[mask_mean_calcite], scatter_data$vr_down[mask_mean_calcite], 
     xlab = "", ylab = "", ylim = c(0.01,100), panel.first = grid(), log = "y",
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
points(scatter_data$latitude[mask_mean_aragonite], scatter_data$vr_down[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
abline(h=1)
lm.out <- lm(log10(scatter_data$vr_down) ~ scatter_data$latitude)
abline(lm.out, lwd = 4, col = "#B2182B")
mtext(text = "latitude  (°N)",side = 1,line = 2)
mtext(text = TeX("$var(\\delta^{18}O_{dw})/var(\\delta^{18}O)$"),side = 2,line = 2)
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
#text(10, 6.7, "calcite", cex = 1.5)
#text(-5, -7., "aragonite", col = "blue", cex = 1.5)
mtext(text = "(c)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#temperature
plot(scatter_data$mean_temp[mask_mean_calcite], scatter_data$vr_down[mask_mean_calcite], 
     yaxt = 'n', xlab = "", ylab = "", ylim = c(0.01,100), log = "y", panel.first = grid(equilogs = FALSE),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
abline(h=1)
points(scatter_data$mean_temp[mask_mean_aragonite], scatter_data$vr_down[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
lm.out <- lm(log(scatter_data$vr_down) ~ scatter_data$mean_temp)
abline(lm.out, lwd = 4, col = "#B2182B")
mtext(text = "Mean temperature (°C)",side = 1,line = 2)
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
mtext(text = "(d)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#geology
plot(c(0,4), c(0.01,100), type = "n", xlab = "", ylab = "", yaxt = 'n', xaxt = 'n', panel.first =grid(equilogs = FALSE), log = "y")
abline(h=1)
boxplot(as.numeric(scatter_data$vr_down[scatter_data$geology == "limestone"]), add = T, at = 0.5, yaxt ="n", axes = F)
boxplot(as.numeric(scatter_data$vr_down[scatter_data$geology == "dolomite"]), add = T, at = 1.5, yaxt ="n", axes = F)
boxplot(as.numeric(scatter_data$vr_down[scatter_data$geology == "marble"]), add = T, at = 2.5, yaxt ="n", axes = F)
boxplot(as.numeric(scatter_data$vr_down[scatter_data$geology == "unknown"]), add = T, at = 3.5, yaxt ="n", axes = F)
mtext(text = "geology", side = 1, line = 2)
axis(1,at=c(0.5,1.5,2.5,3.5),labels=c("limest.", "dolomite", "marble", "?"), cex.axis = 1.1)
mtext(text = "(e)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#prec
plot(scatter_data$mean_prec[mask_mean_calcite], scatter_data$vr_down[mask_mean_calcite], 
     xlab = "", ylab = "", log = "xy", ylim = c(0.01,100), xlim = c(0.2,15), panel.first = grid(equilogs = FALSE),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
abline(h=1)
points(scatter_data$mean_prec[mask_mean_aragonite], scatter_data$vr_down[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
lm.out <- lm(log10(scatter_data$vr_down) ~ scatter_data$mean_prec)
abline(lm.out, lwd = 4, col = "#B2182B")
mtext(text = "Mean annual precipitation (mm/y)",side = 1,line = 2)
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
mtext(text = TeX("$var(\\delta^{18}O_{dw})/var(\\delta^{18}O)$"),side = 2,line = 2)
mtext(text = "(f)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#DJF prec
plot(scatter_data$winter_mean_prec[mask_mean_calcite], scatter_data$vr_down[mask_mean_calcite], 
     yaxt = 'n', xlab = "", ylab = "", xlim = c(0.4,15), log = "xy", ylim = c(0.01,100), yaxt = "n", panel.first = grid(equilogs = FALSE),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
abline(h=1)
points(scatter_data$winter_mean_prec[mask_mean_aragonite], scatter_data$vr_down[mask_mean_aragonite],
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
lm.out <- lm(log10(scatter_data$vr_down) ~ (scatter_data$winter_mean_prec))
abline(lm.out, lwd = 4, col = "#B2182B")
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
mtext(text = "Mean DJF precipitation (mm/y)",side = 1,line = 2)
mtext(text = "(g)", side = 3, line = -1.5, adj = 0.02, cex = 1)


#JJA prec
plot(scatter_data$summer_mean_prec[mask_mean_calcite], scatter_data$vr_down[mask_mean_calcite], 
     yaxt = 'n', xlab = "", ylab = "", xlim = c(0.4,15), log = "xy", ylim = c(0.01,100), yaxt = "n", panel.first = grid(equilogs = FALSE),
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
abline(h=1)
points(scatter_data$summer_mean_prec[mask_mean_aragonite], scatter_data$vr_down[mask_mean_aragonite],
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
lm.out <- lm(log10(scatter_data$vr_down) ~ scatter_data$summer_mean_prec)
abline(lm.out, lwd = 4, col = "#B2182B")
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
mtext(text = "Mean JJA precipitation (mm/y)",side = 1,line = 2)
mtext(text = "(h)", side = 3, line = -1.5, adj = 0.02, cex = 1)

# elevation
plot(scatter_data$elevation[mask_mean_calcite], scatter_data$vr_down[mask_mean_calcite], 
     xlab = "", ylab = "", ylim = c(0.01,100), xlim = c(0,4000), panel.first = grid(), log = "y",
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
abline(h=1)
points(scatter_data$elevation[mask_mean_aragonite], scatter_data$vr_down[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
lm.out <- lm(log10(scatter_data$vr_down) ~ scatter_data$elevation)
abline(lm.out, lwd = 4, col = "#B2182B")
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
mtext(text = "elevation (m)",side = 1,line = 2)
mtext(text = TeX("$var(\\delta^{18}O_{dw})/var(\\delta^{18}O)$"),side = 2,line = 2)
mtext(text = "(i)", side = 3, line = -1.5, adj = 0.02, cex = 1)

plot(scatter_data$elevation_diff[mask_mean_calcite], scatter_data$vr_down[mask_mean_calcite], 
     yaxt = 'n', xlab = "", ylab = "", ylim = c(0.01,100), panel.first = grid(), log = "y",
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
abline(h=1)
points(scatter_data$elevation_diff[mask_mean_aragonite], scatter_data$vr_down[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
lm.out <- lm(log10(scatter_data$vr_down) ~ scatter_data$elevation_diff)
abline(lm.out, lwd = 4, col = "#B2182B")
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
mtext(text = "elevation difference (sim-proxy) (m)",side = 1,line = 2)
mtext(text = "(j)", side = 3, line = -1.5, adj = 0.02, cex = 1)

#cover thickness
plot(scatter_data$cover_thickness[mask_mean_calcite], scatter_data$vr_down[mask_mean_calcite], 
     yaxt = "n", xlab = "", ylab = "", ylim = c(0.01,100), panel.first = grid(), log = "y",
     pch = 16, col = adjustcolor("black", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
abline(h=1)
points(scatter_data$cover_thickness[mask_mean_aragonite], scatter_data$vr_down[mask_mean_aragonite], 
       pch = 16, col = adjustcolor("blue", alpha.f = 0.5), cex = 2, cex.axis = cexaxis)
lm.out <- lm(log10(scatter_data$vr_down) ~ scatter_data$cover_thickness)
abline(lm.out, lwd = 4, col = "#B2182B")
mtext(text = paste0("r²= ", round(summary(lm.out)$r.squared, digits = No.digits)), cex = 1, col = "#B2182B", adj = 1, side = 1, line = -1)
mtext(text = "cover thickness (m)",side = 1,line = 2)
mtext(text = "(k)", side = 3, line = -1.5, adj = 0.02, cex = 1)


dev.off()

rm(scatter_data, mask_mean_aragonite, mask_mean_calcite, data_rec, data_seasonal_a, data_seasonal_b, data_seasonal_c, data_yearly_a, data_yearly_b,
   data_yearly_c, lm.out, plot, cexaxis, entity, ii, No.digits)
