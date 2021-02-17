#################################################
## 3 Paper Plots ################################

library(plyr)
library(dplyr)
library(tidyverse)
library(maps)

# 0.1: Masks

# Prepare masks for different analyses 
mask_mean = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))
mask_var  = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))
mask_spec = logical(length = length(DATA_past1000$CAVES$entity_info$entity_id))

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  data = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
  if(length(data$d18O_measurement) >= 10){mask_mean[ii] = T}
  if(length(data$d18O_measurement) >= 20){mask_var[ii] = T}
  if(length(data$d18O_measurement) >= 30){mask_spec[ii] = T}
}

mask_mean[DATA_past1000$CAVES$entity_info$mineralogy == "mixed"] = F #7 less
mask_var[DATA_past1000$CAVES$entity_info$mineralogy == "mixed"] = F  #7 less
mask_spec[DATA_past1000$CAVES$entity_info$mineralogy == "mixed"] = F #7 less

rm(entity, data, ii)

# 0.2: Cluster

# distance based-clustering
dist<-fossil::earth.dist(cbind(DATA_past1000$CAVES$entity_info$latitude[mask_spec],DATA_past1000$CAVES$entity_info$longitude[mask_spec]),dist=TRUE)
hc<-hclust(dist)
DATA_past1000$CAVES$cluster_list <- list(entity_id = as.numeric(DATA_past1000$CAVES$entity_info$entity_id[mask_spec]), cluster_id = as.numeric(cutree(hc,k=8)))

# manually sort cluster = 9 in south east asia. This is e_ID: 226, 238, 319, 335, 367, 399, 436, 523
for(entity in c(226, 238, 319, 335, 367, 399, 436, 523)){
  if(entity %in% DATA_past1000$CAVES$cluster_list$entity_id){
    DATA_past1000$CAVES$cluster_list$cluster_id[which(DATA_past1000$CAVES$cluster_list$entity_id == entity)] = 9
  }
}

rm(dist, hc, entity)

# 0.3 Gridbox

## Gridbox list

# sorting entities into gridboxes
DATA_past1000$CAVES$gridbox_list <- list(entity_id = as.numeric(DATA_past1000$CAVES$entity_info$entity_id[mask_spec]), gridbox_id = numeric(sum(mask_spec)))

for(ii in 1:sum(mask_spec)){
  e.long = DATA_past1000$CAVES$entity_info$longitude[mask_spec][ii]
  e.lat = DATA_past1000$CAVES$entity_info$latitude[mask_spec][ii]
  DATA_past1000$CAVES$gridbox_list$gridbox_id[ii]  = (ceiling((-1*e.lat+90)/180*73)-1)*96 + ceiling((e.long+180)/3.75)
}
rm(e.lat, e.long, ii)

#################################################
PLOTTING_VARIABLES <- list()
PLOTTING_VARIABLES$WIDTH = 8.3
PLOTTING_VARIABLES$HEIGHT = 5.5

# 1: Timeseries Plot

#Paper Plot for Time Series Analysis. 

# Here: Temp time serie!
source("Functions/Paper_Plots/Fig1_Timeseries.R")

# 2: SISAL location maps

source("Functions/Paper_Plots/Fig2_SISAL.R")
print("entities per cluster")
print(as.tibble(DATA_past1000$CAVES$cluster_list) %>% group_by(cluster_id) %>% count() %>% filter(n>1))
print(paste0("number of entities: ", length(DATA_past1000$CAVES$entity_info$entity_id)))
noC <- DATA_past1000$CAVES$entity_info %>% group_by(site_id) %>% count()
print(paste0("number of entities: ", dim(noC)[1]))
rm(noC)

# 3: Mean State

# Global and cluster mean offset values with 90% confidence intervals
method = "down" # chose analyse-method: "full" for yearly-resolution or "down" for down-sampled to proxy resolution
# Complete Field needed for SD...
source("Functions/Analytics//Mean_Analytics.R")
rm(method)

source("Functions/Paper_Plots/Fig3_MeanState.R")
# Check up drip water equivalent maps 
source("Functions/Checkup_Plots/Fig3_MeanState_CavesOnly.R")

# 4: Scatter Plot

source("Functions/Paper_Plots/Fig4_ScatterPlot.R")

# 5: Variance
#check if variance is changed by drip water conversion
source("Functions/Analytics/Variance_Ratios.R")

source("Functions/Paper_Plots/Fig5_Variance_Map.R")
source("Functions/Paper_Plots/Fig5_Variance_Histo.R")

# 6: Spectrum

#Analysis before plotting!!!
source("Functions/Analytics/Spectra_Analytics.R")

source("Functions/Paper_Plots/Fig6_Spectrum.R")

# 7: Correlation Map
source("Functions/Analytics/Correlation_Analytics.R")

source("Functions/Paper_Plots/Fig7_Correlation.R")

# 8: Network Plot
source("Functions/Analytics/Network_Analytics.R")

source("Functions/Paper_Plots/Fig8_Network.R")
 
# 9: Network Boxplot

source("Functions/Analytics/Ensemble_Analytics.R") #also SF5
source("Functions/Paper_Plots/Fig9_Network_boxplot.R")

