#################################################
## SF Correlation maps with forcing
#################################################

library(plyr)
library(dplyr)
library(tidyverse)

source("Functions/STACYmap_PMIL.R")
source("Functions/aw_mean.R")

#wie muss volcanic forcing gemittelt werden für yearly values?
volcanic_forcing <- Timeseries$volcanic
volcanic_forcing_yearly <- numeric(1251)
volcanic_forcing_yearly[1] <- as.numeric(volcanic_forcing[1])
for(ii in 2:1250){
  year = seq(from = 2,to = 45001, by = 36)[ii]
  volcanic_forcing_yearly[ii] <- mean(as.numeric(volcanic_forcing[year:(year+35)]), na.rm = T)  
}

#plot(volcanic_forcing, xlim = c(-1000,-900))
#points(zoo::zoo(x=volcanic_forcing_yearly,order.by = seq(-1150, 100, length.out = 1251)), col = "blue")

#same time period...
volcanic_forcing_yearly <- volcanic_forcing_yearly[50:1049]
solar_forcing <- Timeseries$solar[1:1000]

#Now calculate correlation maps:
solar_corr <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)), ISOT = array(dim = c(96,73)))
solar_p <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)), ISOT = array(dim = c(96,73)))
volcanic_corr <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)), ISOT = array(dim = c(96,73)))
volcanic_p <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)), ISOT = array(dim = c(96,73)))

for(lon in 1:96){
  for(lat in 2:72){
    for (var in c("TEMP", "PREC", "ISOT")){
      CORR_sol = cor.test(solar_forcing, DATA_past1000$SIM_yearly_a[[var]][lon,lat,])
      CORR_vol = cor.test(volcanic_forcing_yearly, DATA_past1000$SIM_yearly_a[[var]][lon,lat,])
      
      solar_corr[[var]][lon,lat] = CORR_sol$estimate[[1]]
      solar_p[[var]][lon,lat] = CORR_sol$p.value
      volcanic_corr[[var]][lon,lat] = CORR_vol$estimate[[1]]
      volcanic_p[[var]][lon,lat] = CORR_vol$p.value
    }
  }
}

solar_ens_corr <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)),ISOT = array(dim = c(96,73)))
solar_ens_p <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)), ISOT = array(dim = c(96,73)))
volcanic_ens_corr <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)), ISOT = array(dim = c(96,73)))
volcanic_ens_p <- list(TEMP = array(dim = c(96,73)), PREC = array(dim = c(96,73)), ISOT = array(dim = c(96,73)))

for(lon in 1:96){
  for(lat in 2:72){
    for(var in c("TEMP", "PREC", "ISOT")){
      CORR_sol = cor.test(solar_forcing, apply(rbind(DATA_past1000$SIM_yearly_a[[var]][lon,lat,], 
                                                     DATA_past1000$SIM_yearly_b[[var]][lon,lat,],
                                                     DATA_past1000$SIM_yearly_c[[var]][lon,lat,]), 2, mean, na.rm = T))
      CORR_vol = cor.test(volcanic_forcing_yearly, apply(rbind(DATA_past1000$SIM_yearly_a[[var]][lon,lat,], 
                                                               DATA_past1000$SIM_yearly_b[[var]][lon,lat,],
                                                               DATA_past1000$SIM_yearly_c[[var]][lon,lat,]), 2, mean, na.rm = T))
      
      solar_ens_corr[[var]][lon,lat] = CORR_sol$estimate[[1]]
      solar_ens_p[[var]][lon,lat] = CORR_sol$p.value
      volcanic_ens_corr[[var]][lon,lat] = CORR_vol$estimate[[1]]
      volcanic_ens_p[[var]][lon,lat] = CORR_vol$p.value
    }
  }
}

##try maps
#vol_lyr = volcanic_ens_corr
#vol_lyr[volcanic_ens_p>0.1] = NA
#STACYmap(gridlyr = rbind(vol_lyr[49:96,], vol_lyr[1:48,]), centercolor = 0)

## 10% maps
#max_volc_map = apply(DATA_past1000$SIM_yearly_a$ISOT[,,c(which(volcanic_forcing_yearly>0), which(volcanic_forcing_yearly>0)+1)], c(1,2), mean, na.rm = T)
#simpleawmean(max_volc_map - mean_volc_map)

#max_volc_map_TEMP = apply(DATA_past1000$SIM_yearly_a$TEMP[,,which(volcanic_forcing_yearly>0)], c(1,2), mean, na.rm = T)
#mean_volc_map_TEMP = apply(DATA_past1000$SIM_yearly_a$TEMP[,,which(volcanic_forcing_yearly == 0)], c(1,2), mean, na.rm = T)
#simpleawmean(max_volc_map_TEMP - mean_volc_map_TEMP)


CorrMap_Forcing <- list(
  "solar_corr_LM1" = solar_corr,
  "solar_p_LM1" = solar_p,
  "volcanic_corr_LM1" = volcanic_corr,
  "volcanic_p_LM1" = volcanic_p,
  "solar_corr_ens.mean" = solar_ens_corr,
  "solar_p_ens.mean" = solar_ens_p,
  "volcanic_corr_ens.mean" = volcanic_ens_corr,
  "volcanic_p_ens.mean" = volcanic_ens_p, 
  "solar_TS" = solar_forcing,
  "volcanic_TS" = volcanic_forcing_yearly
)

save(CorrMap_Forcing, file = "CorMap_Forcing.RData")

## Correlation maps between a,b,c
CORR_ensemble <- list(
  ISOT = list(
    corr_ab = array(dim = c(96,73)),
    corr_ab_p = array(dim = c(96,73)),
    corr_ac = array(dim = c(96,73)),
    corr_ac_p = array(dim = c(96,73)),
    corr_bc = array(dim = c(96,73)),
    corr_bc_p = array(dim = c(96,73))
    ),
  TEMP = list(
    corr_ab = array(dim = c(96,73)),
    corr_ab_p = array(dim = c(96,73)),
    corr_ac = array(dim = c(96,73)),
    corr_ac_p = array(dim = c(96,73)),
    corr_bc = array(dim = c(96,73)),
    corr_bc_p = array(dim = c(96,73))
  ),
  PREC = list(
    corr_ab = array(dim = c(96,73)),
    corr_ab_p = array(dim = c(96,73)),
    corr_ac = array(dim = c(96,73)),
    corr_ac_p = array(dim = c(96,73)),
    corr_bc = array(dim = c(96,73)),
    corr_bc_p = array(dim = c(96,73))
  )
)


for(lon in 1:96){
  for(lat in 2:72){
    for(var in c("TEMP", "PREC", "ISOT")){
      corr = cor.test(DATA_past1000$SIM_yearly_a[[var]][lon,lat,], DATA_past1000$SIM_yearly_b[[var]][lon,lat,])
      CORR_ensemble[[var]]$corr_ab[lon,lat] = corr$estimate[[1]]
      CORR_ensemble[[var]]$corr_ab_p[lon,lat] = corr$p.value
      
      corr = cor.test(DATA_past1000$SIM_yearly_a[[var]][lon,lat,], DATA_past1000$SIM_yearly_c[[var]][lon,lat,])
      CORR_ensemble[[var]]$corr_ac[lon,lat] = corr$estimate[[1]]
      CORR_ensemble[[var]]$corr_ac_p[lon,lat] = corr$p.value
      
      corr = cor.test(DATA_past1000$SIM_yearly_b[[var]][lon,lat,], DATA_past1000$SIM_yearly_c[[var]][lon,lat,])
      CORR_ensemble[[var]]$corr_bc[lon,lat] = corr$estimate[[1]]
      CORR_ensemble[[var]]$corr_bc_p[lon,lat] = corr$p.value
    }
  }
}

save(CORR_ensemble, file = "CorMap_ensemblesLM1-3.RData")
