#################################################
## Seasonality ensemble increase

library(plyr)
library(dplyr)
library(tidyverse)
library(nest)
library(PaleoSpec)
library(latex2exp)

source("Functions/SubsampleTimeseriesBlock_highresNA.R")

c_season <- list(TEMP = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])),
                 PREC = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])))
p_season <- list(TEMP = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])),
                 PREC = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])))
what_season <- list(TEMP = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])),
                    PREC = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])))
c_yearly <- list(TEMP = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])),
                 PREC = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])))
p_yearly <- list(TEMP = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])),
                 PREC = numeric(length= 3*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])))

for(run in c("a", "b", "c")){
  for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])){
    entity = DATA_past1000$CAVES$entity_info$entity_id[mask_spec][ii]
    data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
    if(run == "a"){jj = ii}
    if(run == "b"){jj = ii + length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])}
    if(run == "c"){jj = ii + 2*length(DATA_past1000$CAVES$entity_info$entity_id[mask_spec])}
    for(var in c("TEMP", "PREC")){
      c_season[[var]][jj] = NA
      p_season[[var]][jj] = NA
      what_season[[var]][jj] = NA
      for(season in c("SPRING", "SUMMER", "AUTUMN", "WINTER")){
        data_season = DATA_past1000$CAVES$season_res[[run]][[season]] %>% filter(entity_id == entity)
        data_sim = SubsampleTimeseriesBlock_highresNA(ts(data = rev(data_season[[var]]), start = LastElement(data_season$year_BP), end = FirstElement(data_season$year_BP)),
                                                      data_rec$interp_age)
        CORR = cor.test(data_sim, data_rec$d18O_measurement)
        if(is.na(p_season[[var]][jj]) & CORR$p.value<0.1){
          c_season[[var]][jj] = CORR$estimate[[1]]
          p_season[[var]][jj] = CORR$p.value
          what_season[[var]][jj] = season
        }else if(abs(CORR$estimate[[1]])>abs(c_season[[var]][jj]) & CORR$p.value < 0.1){
          c_season[[var]][jj] = CORR$estimate[[1]]
          p_season[[var]][jj] = CORR$p.value
          what_season[[var]][jj] = season
        }
      }
      CORR = cor.test(data_rec[[paste0("PREC_", run)]], data_rec$d18O_measurement)
      c_yearly[[var]][jj] = CORR$estimate[[1]]
      p_yearly[[var]][jj] = CORR$p.value
    }
  }
}


#################################################
pdf(file = paste0("Sup_Plots/SF5_SeasonTuning_histo.pdf"), width = 1.3*6, height = 2*PLOTTING_VARIABLES$HEIGHT/1.5)
par(mfrow=c(2,1),oma = c(1,3,0,0) + 0.1,mar = c(3,1,0,1) + 0.1)
hist(c_yearly$TEMP[p_yearly$TEMP<0.1], 
     breaks = 9, border = "white", prob = TRUE, 
     ylim = c(0,3.5), xlim = c(-1,1), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1), cex.axis = 1.5)
lines(density(c_yearly$TEMP[p_yearly$TEMP<0.1], na.rm = T),
      lwd = 2, col = "black")
lines(c(median(c_yearly$TEMP[p_yearly$TEMP<0.1], na.rm = T),
        median(c_yearly$TEMP[p_yearly$TEMP<0.1], na.rm = T)),
      c(0, max(density(c_yearly$TEMP[p_yearly$TEMP<0.1], na.rm = T)$y, na.rm = T)-0.02),
      lwd = 2, col = "black", lty = 2)
lines(density(c_season$TEMP[p_season$TEMP<0.1], na.rm = T),
      lwd = 2, col = "#B2182B")
lines(c(median(c_season$TEMP[p_season$TEMP<0.1], na.rm = T),
        median(c_season$TEMP[p_season$TEMP<0.1], na.rm = T)),
      c(0, max(density(c_season$TEMP[p_season$TEMP<0.1], na.rm = T)$y, na.rm = T)-0.02),
      lwd = 2, col = "#B2182B", lty = 2)

mtext(text = "density",side = 2,line = 2.5, cex = 1.5)
text(-0.5, 2.5, "Season-Tuning", col = "#B2182B", cex = 1.5)
text(-0.5, 2.0, "yearly corr.", col = "black", cex = 1.5)
mtext(text = TeX("$\\rho (T, \\delta^{18}O)$"), side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = 1)
text(1, 2.75, paste0("DJF: ", floor(sum(what_season$TEMP == "WINTER", na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
text(1, 2.25, paste0("MAM: ", floor(sum(what_season$TEMP == "SPRING", na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
text(1, 1.75, paste0("JJA: ", floor(sum(what_season$TEMP == "SUMMER", na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
text(1, 1.25, paste0("SON: ", floor(sum(what_season$TEMP == "AUTUMN", na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
text(1, 0.75, paste0("Year: ", floor(sum(p_yearly$TEMP<0.1, na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
mtext(text = "(a)", side = 3, line = -2, adj = 0, col = "black", cex = 1.5, at =-1)


hist(c_yearly$PREC[p_yearly$PREC<0.1],
     breaks = 9, border = "white", prob = TRUE,
     ylim = c(0,3.5), xlim = c(-1,1), xlab = "",xaxt = 'n', ylab = "",
     main = "", cex.main = 1.5, cex.axis = 1.5)
axis(side = 1, at = c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1), cex.axis = 1.5)
lines(density(c_yearly$PREC[p_yearly$PREC<0.1], na.rm = T),
      lwd = 2, col = "black")
lines(c(median(c_yearly$PREC[p_yearly$PREC<0.1], na.rm = T),
        median(c_yearly$PREC[p_yearly$PREC<0.1], na.rm = T)),
      c(0, max(density(c_yearly$PREC[p_yearly$PREC<0.1], na.rm = T)$y, na.rm = T)-0.02),
      lwd = 2, col = "black", lty = 2)
lines(density(c_season$PREC[p_season$PREC<0.1], na.rm = T),
      lwd = 2, col = "#B2182B")
lines(c(median(c_season$PREC[p_season$PREC<0.1], na.rm = T),
        median(c_season$PREC[p_season$PREC<0.1], na.rm = T)),
      c(0, max(density(c_season$PREC[p_season$PREC<0.1], na.rm = T)$y, na.rm = T)-0.02),
      lwd = 2, col = "#B2182B", lty = 2)
mtext(text = TeX("$\\rho (P, \\delta^{18}O)$"), side = 3, line = -2, adj = 1,col = "black", cex = 1.5, at = 1)
text(1, 2.75, paste0("DJF: ", floor(sum(what_season$PREC == "WINTER", na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
text(1, 2.25, paste0("MAM: ", floor(sum(what_season$PREC == "SPRING", na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
text(1, 1.75, paste0("JJA: ", floor(sum(what_season$PREC == "SUMMER", na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
text(1, 1.25, paste0("SON: ", floor(sum(what_season$PREC == "AUTUMN", na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
text(1, 0.75, paste0("Year: ", floor(sum(p_yearly$PREC<0.1, na.rm = T)/3)), col = "black", cex = 1.5, adj = 1)
mtext(text = "(b)", side = 3, line = -2, adj = 0, col = "black", cex = 1.5, at =-1)

dev.off()

rm(c_season, c_yearly, CORR, data_rec, data_season, p_season, p_yearly, what_season, data_sim, entity, ii, jj, run, season, var)
