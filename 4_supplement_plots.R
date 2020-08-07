#################################################
## 4 Supplement plots ###########################

library(plyr)
library(dplyr)
library(tidyverse)
library(maps)
library(nest)
library(PaleoSpec)

# SF1: All Time Series with offsets of records

source("Functions/Supplement_Plots/SF_All_Timeseries.R")

# SF43: Vegetation

source("Functions/Supplement_Plots/Vegetation_Plot.R")

# SF4: Variance Scatter Plot

source("Functions/Supplement_Plots/SF_ScatterPlot_Variance.R")

# SF5: All Spectra of records. 

source("Functions/Supplement_Plots/SF15_AllSpectra.R")

# SF7: Network Plot with yearly resolution
#--> already in paper plot 8!!!
print("In Functions/Paper_Plots/Fig8_Network.R")

# SF8: Network maps and matrix comparison DS Europe, China, South+North America

source("Functions/Supplement_Plots/SF_Network_ClusterMap.R")

# SF12: Time scale dependent variance histogram

source("Functions/Supplement_Plots/SF2_TimeScaleVariance.R")

# SF13: Season correlation map with strongest season
# SF15: Season Tuned correlation map

source("Functions/Supplement_Plots/SF_Seasonality_Map1.R")#<- both maps

# SF14: Correlation Tuning Change compare to before

source("Functions/Analytics/Season_CorrTuning.R") #also SF5

# SF16: Network tuning via closest cave and smallest offset

source("Functions/Supplement_Plots/SF_Network_Tuning.R")







