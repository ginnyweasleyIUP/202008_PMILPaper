#################################################
## DISCUSSION Time Scale Variance ###############
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(PaleoSpec)
library(nest)

source('Functions/Filter/EASY_Sensor_WM4.R')
source('Functions/Filter/filter_function3.R')
source('Functions/filter_window.R')

VARIANCE <- list(rec_short = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 rec_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 rec_tot = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 rec_tot_spec = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_full_short = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_full_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_full_tot = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_down_short = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_down_long = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_down_tot = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_down_long_spec = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_down_short_spec = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_down_tot_spec = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_down_filter_long_spec = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_down_filter_short_spec = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])),
                 sim_down_filter_tot_spec = numeric(length = length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])))

for(entity_number in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_var])){
  TS = list()
  entity = DATA_past1000$CAVES$entity_info$entity_id[mask_var][entity_number]
  if(entity %in% c(115, 124, 128, 163, 319, 405, 523)){
    for(list in 1:length(VARIANCE)){
      VARIANCE[[list]][entity_number] = NA
    }
    next
  }
  
  data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
  data_yearly = DATA_past1000$CAVES$yearly_res$a %>% filter(entity_id == entity)
  diff = floor((LastElement(data_rec$interp_age)-FirstElement(data_rec$interp_age))/length(data_rec$interp_age))
  if(diff<1){diff = 1}
  
  TS$TS_rec = zoo(x = data_rec$d18O_dw_eq_a, order.by = data_rec$interp_age)
  TS$TS_rec_eq = PaleoSpec::MakeEquidistant(data_rec$interp_age, data_rec$d18O_measurement,
                                            time.target = seq(from = FirstElement(data_rec$interp_age), 
                                                              to = LastElement(data_rec$interp_age), by = diff))
  TS$TS_sim_full = zoo(x = data_yearly$ITPC, order.by = data_yearly$year_BP)
  TS$TS_sim_down = zoo(x = data_rec$ITPC_a, order.by = data_rec$interp_age)
  TS$TS_sim_down_eq = PaleoSpec::MakeEquidistant(data_rec$interp_age, data_rec$ITPC_a,
                                                 time.target = seq(from = FirstElement(data_rec$interp_age), 
                                                                   to = LastElement(data_rec$interp_age), by = diff))
  
  Results <- easy_sensor_wm4(diff, na.omit(rev(TS$TS_sim_down_eq)), 3.0)
  time_new = seq(LastElement(data_rec$interp_age), LastElement(data_rec$interp_age)-diff*length(Results), by = -1*diff)
  TS$TS_sim_down_filter = filter_window(ts(data = rev(Results), start = LastElement(time_new), end = FirstElement(time_new), deltat = diff))
  
  for(method in c("rec", "sim_full", "sim_down")){
    for(tsc_length in c("long", "short")){
      if(tsc_length == "long"){
        tsc = c(50,200)
      }else{
        tsc = c(10,50)
      }
      if(is.na(tsc_dep_var(timser = TS[[paste0("TS_", method)]], tsc.in = tsc)$var.tsc)){
        VARIANCE[[paste0(method, "_", tsc_length)]][entity_number] = NA
      }else{
        temp = tsc_dep_var(timser = TS[[paste0("TS_", method)]], tsc.in = tsc)
        VARIANCE[[paste0(method, "_", tsc_length)]][entity_number]  = temp$var.tsc/temp$var.tot
      }
    }
    VARIANCE[[paste0(method,"_tot")]][entity_number] = var(TS[[paste0("TS_", method)]], na.rm = T)
  }
  
  SPEC<- list()
  SPEC$Spec_rec = spectrum(na.omit(TS$TS_rec_eq), plot = F)
  SPEC$Spec_sim_down = spectrum(as.ts(na.omit(TS$TS_sim_down_eq)), plot = F)
  SPEC$Spec_sim_down_filter  = spectrum(na.omit(TS$TS_sim_down_filter), plot = F)
  
  for(method in c("sim_down", "sim_down_filter")){
    #short(<50y)
    if(LastElement(SPEC[[paste0("Spec_",method)]]$freq)<=1/50){
      VARIANCE[[paste0(method, "_short_spec")]][entity_number] = NA
    }else{
      VARIANCE[[paste0(method, "_short_spec")]][entity_number] = GetVarFromSpectra(SPEC[[paste0("Spec_",method)]],
                                                                                                      f= c(1/50, LastElement(SPEC[[paste0("Spec_",method)]]$freq)))$var
    }
    #long(>50y)
    if(FirstElement(SPEC[[paste0("Spec_",method)]]$freq)>=1/50){
      VARIANCE[[paste0(method, "_long_spec")]][entity_number] = NA
    }else{
      VARIANCE[[paste0(method, "_long_spec")]][entity_number] = GetVarFromSpectra(SPEC[[paste0("Spec_",method)]],
                                                                                                      f= c(FirstElement(SPEC[[paste0("Spec_",method)]]$freq),
                                                                                                           1/50))$var
    }
    #total
    VARIANCE[[paste0(method, "_tot_spec")]][entity_number] = GetVarFromSpectra(SPEC[[paste0("Spec_",method)]],
                                                                               f= c(FirstElement(SPEC[[paste0("Spec_",method)]]$freq),
                                                                                    LastElement(SPEC[[paste0("Spec_",method)]]$freq)))$var
  }
  
  VARIANCE$rec_tot_spec[entity_number] = GetVarFromSpectra(SPEC$Spec_rec, f= c(FirstElement(SPEC$Spec_rec$freq),LastElement(SPEC$Spec_rec$freq)))$var
}

  

cex_text = 1.5
cex_axis = 2

for(plot in 1:1){
  pdf(file = paste0("Sup_Plots/SF2_Variance_TimeScaleDep.pdf"), width = 1.3*6, height = 3*PLOTTING_VARIABLES$HEIGHT/1.5)
  par(mfrow=c(3,1),oma = c(1,3,0,0) + 0.1,mar = c(3,1,0,1) + 0.1)
  hist(log10(VARIANCE$rec_short), 
       breaks = 9, border = "white", prob = TRUE, 
       ylim = c(0,6), xlim = c(-2,0), xlab = "",xaxt = 'n', ylab = "",
       main = "", cex.main = cex_text, cex.axis = cex_axis)
  axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0), 
       labels = c(0.001,0.01, 0.1, 1), cex.axis = cex_axis)
  lines(density(log10(VARIANCE$rec_short), na.rm = T),
        lwd = 2, col = "black")
  lines(c(median(log10(VARIANCE$rec_short), na.rm = T),median(log10(VARIANCE$rec_short), na.rm = T)), 
        c(0, max(density(log10(VARIANCE$rec_short), na.rm = T)$y)-0.01),
        lwd = 2, col = "black", lty = 2)
  lines(density(log10(VARIANCE$rec_long), na.rm = T),
        lwd = 2, col = "#B2182B")
  lines(c(median(log10(VARIANCE$rec_long), na.rm = T),median(log10(VARIANCE$rec_long), na.rm = T)), 
        c(0, max(density(log10(VARIANCE$rec_long), na.rm = T)$y)-0.01),
        lwd = 2, col = "#B2182B", lty = 2)
  mtext(text = "density",side = 2,line = 2.5, cex = cex_text)
  text(log10(0.1), 4, "long (50-200y)", col = "#B2182B", cex = cex_axis)
  text(log10(0.1), 3, "short (10-50y)", col = "black", cex = cex_axis)
  
  text(log10(0.01), 4, paste0("#Records: ", length(na.omit(VARIANCE$rec_long))), col = "#B2182B", cex = cex_axis, adj = 0)
  text(log10(0.01), 3, paste0("#Records: ", length(na.omit(VARIANCE$rec_short))), col = "black", cex = cex_axis, adj = 0)
  
  mtext(text = "Record", side = 3, line = -2, adj = 1,col = "black", cex = cex_text, at = log10(1))
  mtext(text = "(a)", side = 3, line = -2, adj = 0,col = "black", cex = cex_text, at = log10(0.01))
  
  hist(log10(VARIANCE$sim_full_short), 
       breaks = 9, border = "white", prob = TRUE, 
       ylim = c(0,6), xlim = c(-2,0), xlab = "",xaxt = 'n', ylab = "",
       main = "", cex.main = cex_text, cex.axis = cex_axis)
  axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0), 
       labels = c(0.001,0.01, 0.1, 1), cex.axis = cex_axis)
  lines(density(log10(VARIANCE$sim_full_short), na.rm = T),
        lwd = 2, col = "black")
  lines(c(median(log10(VARIANCE$sim_full_short), na.rm = T),median(log10(VARIANCE$sim_full_short), na.rm = T)), 
        c(0, max(density(log10(VARIANCE$sim_full_short), na.rm = T)$y)-0.01),
        lwd = 2, col = "black", lty = 2)
  lines(density(log10(VARIANCE$sim_full_long), na.rm = T),
        lwd = 2, col = "#B2182B")
  lines(c(median(log10(VARIANCE$sim_full_long), na.rm = T),median(log10(VARIANCE$sim_full_long), na.rm = T)), 
        c(0, max(density(log10(VARIANCE$sim_full_long), na.rm = T)$y)-0.01),
        lwd = 2, col = "#B2182B", lty = 2)
  mtext(text = "density",side = 2,line = 2.5, cex = cex_text)
  mtext(text = "Sim full", side = 3, line = -2, adj = 1,col = "black", cex = cex_text, at = log10(1))
  mtext(text = "(b)", side = 3, line = -2, adj = 0,col = "black", cex = cex_text, at = log10(0.01))
  text(log10(0.01), 4, paste0("#Records: ", length(na.omit(VARIANCE$sim_full_long))), col = "#B2182B", cex = cex_axis, adj = 0)
  text(log10(0.01), 3, paste0("#Records: ", length(na.omit(VARIANCE$sim_full_short))), col = "black", cex = cex_axis, adj = 0)
  
  hist(log10(VARIANCE$sim_down_short), 
       breaks = 9, border = "white", prob = TRUE, 
       ylim = c(0,6), xlim = c(-2,0), xlab = "",xaxt = 'n', ylab = "",
       main = "", cex.main = cex_text, cex.axis = cex_axis)
  axis(side = 1, at = c(log10(0.001), log10(0.01), log10(0.1), 0), 
       labels = c(0.001,0.01, 0.1, 1), cex.axis = cex_axis)
  lines(density(log10(VARIANCE$sim_down_short), na.rm = T),
        lwd = 2, col = "black")
  lines(c(median(log10(VARIANCE$sim_down_short), na.rm = T),median(log10(VARIANCE$sim_down_short), na.rm = T)), 
        c(0, max(density(log10(VARIANCE$sim_down_short), na.rm = T)$y)-0.01),
        lwd = 2, col = "black", lty = 2)
  lines(density(log10(VARIANCE$sim_down_long), na.rm = T),
        lwd = 2, col = "#B2182B")
  lines(c(median(log10(VARIANCE$sim_down_long), na.rm = T),median(log10(VARIANCE$sim_down_long), na.rm = T)), 
        c(0, max(density(log10(VARIANCE$sim_down_long), na.rm = T)$y)-0.01),
        lwd = 2, col = "#B2182B", lty = 2)
  mtext(text = "density",side = 2,line = 2.5, cex = cex_text)
  mtext(text = "Sim down", side = 3, line = -2, adj = 1,col = "black", cex = cex_text, at = log10(1))
  mtext(text = "var.tsc(d18Opw)/var.tot(d18Opw)",side = 1,line = 2.7, cex = cex_text)
  mtext(text = "(c)", side = 3, line = -2, adj = 0,col = "black", cex = cex_text, at = log10(0.01))
  
  text(log10(0.01), 4, paste0("#Records: ", length(na.omit(VARIANCE$sim_down_long))), col = "#B2182B", cex = cex_axis, adj = 0)
  text(log10(0.01), 3, paste0("#Records: ", length(na.omit(VARIANCE$sim_down_short))), col = "black", cex = cex_axis, adj = 0)
  
  
  dev.off()
  
  
}
for(analysis in 1:1){
  # Analytics we need for the paper:
  No.digits = 2
  
  # the median resolution of the records of $10.1 \unit{y} (8.39\unit{y}, 11.73\unit{y})$
  res = numeric(length(DATA_past1000$CAVES$entity_info$entity_id[mask_var]))
  for(ii in 1:length(res)){
    entity = DATA_past1000$CAVES$entity_info$entity_id[mask_var][ii]
    data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
    res[ii] = median(diff(data_rec$interp_age), na.rm = T)
  }
  bstrap <- c()
  for (i in 1:1000){
    bstrap <- c(bstrap,median(sample(res,length(res),replace=T), na.rm = T))
  }
  
  print(paste0("the median resolution of the records is $", round(median(res, na.rm = T), digits = No.digits), " (", round(quantile(bstrap,0.05), digits = No.digits), ", ", round(quantile(bstrap,0.95), digits = No.digits), ")$"))
  
  # On shorter timescales ($<50y$), the filter reduces the variance to 0.41 (0.21, 7.24)} of the original variance
  bstrap <- c()
  for (i in 1:1000){
    bstrap <- c(bstrap,median(sample(VARIANCE$sim_down_filter_short_spec/VARIANCE$sim_down_short_spec,length(VARIANCE$sim_down_short_spec),replace=T), na.rm = T))
  }
  print(paste0("On shorter timescales(<50y), the filter reduces the variance to ", round(median(VARIANCE$sim_down_filter_short_spec/VARIANCE$sim_down_short_spec, na.rm = T), digits = No.digits),
               " (", round(quantile(bstrap,0.05), digits = No.digits), ", ", round(quantile(bstrap,0.95), digits = No.digits), ") of the unfiltered variance."))
  
  # while it increases the variance to more than 5.66 (1.11, 150.00)} of the original down-sampled variance on longer time scales ($>50y$).
  bstrap <- c()
  for (i in 1:1000){
    bstrap <- c(bstrap,median(sample(VARIANCE$sim_down_filter_long_spec/VARIANCE$sim_down_long_spec,length(VARIANCE$sim_down_long_spec),replace=T), na.rm = T))
  }
  print(paste0("On longer timescales(<50y), the filter reduces the variance to ", round(median(VARIANCE$sim_down_filter_long_spec/VARIANCE$sim_down_long_spec, na.rm = T), digits = No.digits),
               " (", round(quantile(bstrap,0.05), digits = No.digits), ", ", round(quantile(bstrap,0.95), digits = No.digits), ") of the unfiltered variance."))
  
  
  
  #The total variance of the karst-filtered down-sampled simulated \deltaO increases to $3.15\times (0.54, 150.00)$} the unfiltered down-sampled varianc
  bstrap <- c()
  for (i in 1:1000){
    bstrap <- c(bstrap,median(sample(VARIANCE$sim_down_filter_tot_spec/VARIANCE$sim_down_tot_spec,length(VARIANCE$sim_down_tot_spec),replace=T), na.rm = T))
  }
  print(paste0("The filter reduces the total variance to ", round(median(VARIANCE$sim_down_filter_tot_spec/VARIANCE$sim_down_tot_spec, na.rm = T), digits = No.digits),
               " (", round(quantile(bstrap,0.05), digits = No.digits), ", ", round(quantile(bstrap,0.95), digits = No.digits), ") of the unfiltered variance."))
  
  
  #which equals $2.23\times(0.27, 35.70)$} the record variance.
  bstrap <- c()
  for (i in 1:1000){
    bstrap <- c(bstrap,median(sample(VARIANCE$sim_down_filter_tot_spec/VARIANCE$rec_tot_spec,length(VARIANCE$rec_tot_spec),replace=T), na.rm = T))
  }
  print(paste0("The filtered total varians equals ", round(median(VARIANCE$sim_down_filter_tot_spec/VARIANCE$rec_tot_spec, na.rm = T), digits = No.digits),
               " (", round(quantile(bstrap,0.05), digits = No.digits), ", ", round(quantile(bstrap,0.95), digits = No.digits), ") of the total record variance."))
  
  
  
  #The ratio between record variance and down-sampled variance compared to the median in Fig. \ref{fig:4_Variance} is $2.26\times(0.26, 22.81)$ 
  # which is 25\% higher. 
  bstrap <- c()
  for (i in 1:1000){
    bstrap <- c(bstrap,median(sample(VARIANCE$rec_tot_spec/VARIANCE$sim_down_tot_spec,length(VARIANCE$rec_tot_spec),replace=T), na.rm = T))
  }
  bstrap_2 <- c()
  for (i in 1:1000){
    bstrap_2 <- c(bstrap_2,median(sample(((VARIANCE$rec_tot_spec/VARIANCE$sim_down_tot_spec)/1.97-1)*100,length(VARIANCE$rec_tot_spec),replace=T), na.rm = T))
  }
  print(paste0("The total record variance equals ", round(median(VARIANCE$rec_tot_spec/VARIANCE$sim_down_tot_spec, na.rm = T), digits = No.digits),
               " (", round(quantile(bstrap,0.05), digits = No.digits), ", ", round(quantile(bstrap,0.95), digits = No.digits), ") of the total downsampled variance,",
               " which is ", round((median(VARIANCE$rec_tot_spec/VARIANCE$sim_down_tot_spec, na.rm = T)/1.97-1)*100, digits = No.digits),
               " (", round(quantile(bstrap_2,0.05), digits = No.digits), ", ", round(quantile(bstrap_2,0.95), digits = No.digits), ") % larger than in the variance analysis." ))
  
  
}

rm(data_rec, data_yearly, entities, SPEC, temp, TS, analysis, bstrap, bstrap_2, cex_axis, cex_text, diff, entity,
   entity_number, GAUSS, i, ii, list, method, No.digits, plot, res, Results, time_new, tsc, tsc_length)
