## Spectrum Plot 

## PLOT #########################################

filter <- list(
  ISOT = c(2.0, 16.0, 8.0),
  ITPC = c(3.0, 9.0, 3.0)
)
cex_text = 1.0

COLZ <- c("#1A2254", "#0A4296", "#278BCE","#91002B", "#BD6B73", "black")
COLZ <- c("#1A2254", "#0A4296", "#278BCE","#91002B", "#D6604D", "black")


#################################################
## WITH FILTER ##################################
#################################################

source('Functions/Filter/EASY_Sensor_WM4.R')
source('Functions/Filter/filter_function3.R')

entity = 240
data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
data_yearly = DATA_past1000$CAVES$yearly_res$a %>% filter(entity_id == entity)

for(var in c("ISOT", "ITPC")){
  #Filter 
  #full -> down
  
  Results <- easy_sensor_wm4(1.0, na.omit(data_yearly[[var]]), filter[[var]][1])#filter[[var]][1])
  time_new = seq(FirstElement(data_yearly$year_BP), LastElement(data_yearly$year_BP)-(length(Results)-length(data_yearly$ITPC)), by = -1)
  TS_full_down  = ts(data = rev(Results), start = LastElement(time_new), end = FirstElement(time_new), deltat = 1)

  #full-> rec
  Results <- easy_sensor_wm4(1.0, na.omit(data_yearly[[var]]), filter[[var]][2])
  time_new = seq(FirstElement(data_yearly$year_BP), LastElement(data_yearly$year_BP)-(length(Results)-length(data_yearly$ITPC)), by = -1)
  TS_full_rec  = ts(data = rev(Results), start = LastElement(time_new), end = FirstElement(time_new), deltat = 1)
  
  #down->rec
  diff = floor((LastElement(data_rec$interp_age)-FirstElement(data_rec$interp_age))/length(data_rec$interp_age))
  if(diff<1){diff = 1}
  record <- as.numeric(PaleoSpec::MakeEquidistant(data_rec$interp_age,
                                       data_rec[[paste0(var,"_a")]],
                                       time.target = seq(from = FirstElement(data_rec$interp_age), to = LastElement(data_rec$interp_age), by = diff)))
  Results <- easy_sensor_wm4(diff, rev(na.omit(record)), filter[[var]][3])
  time_new = seq(LastElement(data_rec$interp_age), LastElement(data_rec$interp_age)-diff*length(Results)+1, by = -1*diff)
  TS_down_rec = ts(data = rev(Results), start = LastElement(time_new), end = FirstElement(time_new), deltat = diff)
  
  cairo_pdf(file = paste0("Paper_Plots/Fig6_Spectra_Compare_Filter_",var,".pdf"), width = PLOTTING_VARIABLES$WIDTH*2/3, height = PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  #par(mfrow=c(1,2))
  cex_text = 1
  cex_axis = 0.6
  cex_axis_text = 0.7
  par(mar=c(1,0,0,1),oma=c(3,5,2,3),xaxs="i",yaxs="i",cex=1,lwd=2)
  layout(matrix(c(1,2,3,4,4,4), 3, 2, byrow = F))
  
  plot(data_rec$interp_age, data_rec$d18O_dw_eq_a, col = "black", type = "l", ylim = c(-10,-5.0), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col=COLZ[6])
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),las=1,col=COLZ[6], cex = cex_axis)
  text(1050, -5.5, TeX("Example: eID 240 (Bunker Cave)"), col = COLZ[6], adj = 1, cex = cex_text)
  mtext("(a)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  ######################################################
  ## der Filter funktioniert nicht! 
  ## die Zeitachse des Filters ist falschrum!
  ## kontrolliere auch alle Filter in der Spectral Analysis...
  
  
  plot(data_yearly$year_BP, data_yearly[[var]], 
       col = COLZ[1], type = "l", ylim = c(-10,-5.0), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  lines(window(TS_full_down, start = 110, end = 1100), col = COLZ[2])
  lines(window(TS_full_rec, start = 100, end = 1050), col = COLZ[3])
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col=COLZ[1])
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),col=COLZ[1], cex = cex_axis)
  text(1050, -5.5, TeX("HadCM3 at eID 240"), col = COLZ[1], adj = 1, cex = cex_text)
  text(1050, -6.0, paste0(filter[[var]][1], "y filter"), col = COLZ[2], adj = 1, cex = cex_text)
  text(1050, -9.5, paste0(filter[[var]][2], "y filter"), col = COLZ[3], adj = 1, cex = cex_text)
  mtext(TeX("$\\delta^{18}O_{pw}$"),side=2, at = -8, line = 2.5, las = 1, col = COLZ[1], cex = cex_axis_text, las = 0)
  mtext(text = "[‰]", side = 2, at = -6.5, line = 2.5, las = 1, col = COLZ[1], cex = cex_axis_text, las = 0)
  mtext("(b)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  plot(data_rec$interp_age, data_rec[[paste0(var,"_a")]],
       col = COLZ[4], type = "l", ylim = c(-10,-5.0), xlim = c(0,1100), xaxt = "n", yaxt = "n",ylab = "", xlab = "")
  lines(window(TS_down_rec, start = 140, end = 1100), col = COLZ[5])
  axis(2,at=seq(-10,-6,by=1),labels=FALSE,col=COLZ[4])
  mtext(side=2,at=seq(-10,-6,by=1),line = 1, seq(-10,-6,by=1),las=1,col=COLZ[4], cex = cex_axis)
  text(1050, -5.5, TeX("HadCM3 at eID 240 down-sampled"), col = COLZ[4], adj = 1, cex = cex_text)
  text(1050, -6.0, paste0(filter[[var]][3], "y filter"), col = COLZ[5], adj = 1, cex = cex_text)
  mtext("(c)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  
  
  axis(1,at=seq(0,1000,by=200),labels=FALSE,col="black")
  mtext(side=1,at=seq(0,1000,by=200),line = 1, seq(0,1000,by=200),las=1,col="black", cex = cex_axis)
  mtext("years BP", side = 1, line = 2, at = 500, cex = cex_axis_text)
  
  #png(file = paste0("Plots/Paper_Plot_7_Spectra_",var,".png"), width = 70*PLOTTING_VARIABLES$WIDTH*2/3, height = 70*PLOTTING_VARIABLES$HEIGHT*1.2*2/3)
  LPlot(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_",var)]]$spec, col = adjustcolor(COLZ[1], alpha = 0.5), 
        ylim = c(0.00005,1000), xlim = c(1/300, 0.5),
        ylab = "",
        xaxt = 'n',
        yaxt = "n",
        xlab = "", lwd = 2)
  LLines(LogSmooth(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_",var)]]$spec), col = COLZ[1], lw = 3)#,
  #main = TeX("Mean Spectra from cave locations (res>8)"))
  mtext("Period (y)", side = 1, line= 2, cex = cex_axis_text)
  mtext(expression("PSD [‰"^"2"*"/year]"), side = 4, line= 2, at = 1, las = 0, cex = cex_axis_text)
  LLines(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_down_",var)]]$spec, col = COLZ[2], lty = 3, lw = 2)
  LLines(LogSmooth(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_down_",var)]]$spec), col = COLZ[2], lw = 2)
  LLines(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_rec_",var)]]$spec, col = COLZ[3], lty = 3, lw = 2)
  LLines(LogSmooth(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_full_rec_",var)]]$spec), col = COLZ[3], lw = 2)
  #text(0.2, 8e2, "HadCM3 yearly res.", col = "#074893")
  #text(0.3, 3e2, "5y filter", col = "#074893")
  #text(0.3, 1e2, "50y filter", col = "#074893")
  
  axis(side = 1, at = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.2, 0.5), labels = FALSE)
  mtext(side=1,at=c(0.005, 0.01, 0.02, 0.05, 0.2, 0.5),line = 1, 
        c(1/0.005, 1/0.01, 1/0.02, 1/0.05, 1/0.2, 1/0.5),las=1,col="black", cex = cex_axis)
  axis(side = 4, at = c(1e-3, 1e-1, 1e1, 1e3), labels = FALSE)
  mtext(side=4,at=c(1e-3, 1e-1, 1e1, 1e3),line = 1, c(1e-3, 0.1, 10, 1000),las=1,col="black", cex = cex_axis)
  
  LLines(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_ds_",var)]]$spec, col = adjustcolor(COLZ[4], alpha = 0.5), lw = 2)
  LLines(LogSmooth(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_ds_",var)]]$spec), col = COLZ[4], lw = 2)
  LLines(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_down_rec_",var)]]$spec, col = COLZ[5], lty = 3, lw = 2)
  LLines(LogSmooth(SPECTRA$MEAN_SPEC_WEIGH[[paste0("SIM_down_rec_",var)]]$spec), col = COLZ[5], lw = 2)
  #text(0.02, 0.01, "HadCM3 down-sampled", col = COLZ[4])
  #text(0.02, 0.001, "10y filter", col = COLZ[4])
  
  LLines(SPECTRA$MEAN_SPEC_WEIGH$Record$spec, col = adjustcolor(COLZ[6], alpha = 0.5), lw = 2)
  LLines(LogSmooth(SPECTRA$MEAN_SPEC_WEIGH$Record$spec), col = COLZ[6], lw = 3)
  #text(0.02, 0.0005, "Records", col = "black")
  
  legend("bottomleft", legend = c("HadCM3 (full)", paste0(" ...  ",filter[[var]][1],"y filter"), paste0(" ...  ",filter[[var]][2],"y filter"), "HadCM3 (down-sampled)", paste0(" ...  ",filter[[var]][3],"y filter"), "Records"), 
         col = c(COLZ[1],COLZ[2],COLZ[3],COLZ[4],COLZ[5],COLZ[6]), lwd = c(2,2,2,2,2,3), lty = c(1,1,1,1,1,1), bty = "n", cex = cex_text)
  
  if(var == "ISOT"){
    text(0.4, 500, TeX("Spectrum HadCM3: $\\delta^{18}O$"), adj = 1, cex = cex_text)
  }
  if(var == "ITPC"){
    text(0.4, 500, TeX("Spectrum HadCM3: $\\delta^{18}O_{pw}$"), adj = 1, cex = cex_text)
  }
  mtext("d)", side = 3, line = -1.3, adj = 0.025, cex = cex_axis_text)
  
  
  dev.off()
  
}

rm(cex_axis, cex_axis_text, cex_text, COLZ, diff, Results, time_new, data_rec, entity, var, record)
rm(TS_down_rec, TS_full_down, TS_full_rec, filter, data_yearly)
