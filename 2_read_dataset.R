#################################################
## 2 Read in Dataset ############################

# Fields that we need:
# TEMP MEAN
# PREC MEAN
# ISOT MEAN
# ITPC MEAN
# COR MAP TEMP-ITPC
# COR MAP PREC-ITPC

DATA_past1000 <- list()
DATA_past1000$CAVES <- list()
DATA_past1000$CAVES$entity_info <- read.csv("Data/SISAL_HadCM3_entity_info.csv")

DATA_past1000$CAVES$entity_info$elevation <- as.numeric(DATA_past1000$CAVES$entity_info$elevation)
DATA_past1000$CAVES$entity_info$geology <- as.character(DATA_past1000$CAVES$entity_info$geology)
DATA_past1000$CAVES$entity_info$mineralogy <- as.character(DATA_past1000$CAVES$entity_info$mineralogy)
DATA_past1000$CAVES$entity_info$distance_entrance <- as.numeric(levels(DATA_past1000$CAVES$entity_info$distance_entrance))[DATA_past1000$CAVES$entity_info$distance_entrance]
DATA_past1000$CAVES$entity_info$cover_thickness <- as.numeric(levels(DATA_past1000$CAVES$entity_info$cover_thickness))[DATA_past1000$CAVES$entity_info$cover_thickness]

DATA_past1000$CAVES$record_res <- read.csv("Data/SISAL_HadCM3_ds.csv")
DATA_past1000$CAVES$yearly_res <- list()
DATA_past1000$CAVES$yearly_res$a <- read.csv("Data/SISAL_HadCM3_xnapa_PMIL_yearly.csv")
DATA_past1000$CAVES$yearly_res$b <- read.csv("Data/SISAL_HadCM3_xnapb_PMIL_yearly.csv")
DATA_past1000$CAVES$yearly_res$c <- read.csv("Data/SISAL_HadCM3_xnapc_PMIL_yearly.csv")
DATA_past1000$CAVES$season_res <- list()
DATA_past1000$CAVES$season_res$a <- list(
  WINTER = read.csv("Data/Seasonal/SISAL_HadCM3_xnapa_PMIL_WINTER.csv"),
  SPRING = read.csv("Data/Seasonal/SISAL_HadCM3_xnapa_PMIL_SPRING.csv"),
  SUMMER = read.csv("Data/Seasonal/SISAL_HadCM3_xnapa_PMIL_SUMMER.csv"),
  AUTUMN = read.csv("Data/Seasonal/SISAL_HadCM3_xnapa_PMIL_AUTUMN.csv")
)
DATA_past1000$CAVES$season_res$b <- list(
  WINTER = read.csv("Data/Seasonal/SISAL_HadCM3_xnapb_PMIL_WINTER.csv"),
  SPRING = read.csv("Data/Seasonal/SISAL_HadCM3_xnapb_PMIL_SPRING.csv"),
  SUMMER = read.csv("Data/Seasonal/SISAL_HadCM3_xnapb_PMIL_SUMMER.csv"),
  AUTUMN = read.csv("Data/Seasonal/SISAL_HadCM3_xnapb_PMIL_AUTUMN.csv")
)
DATA_past1000$CAVES$season_res$c <- list(
  WINTER = read.csv("Data/Seasonal/SISAL_HadCM3_xnapc_PMIL_WINTER.csv"),
  SPRING = read.csv("Data/Seasonal/SISAL_HadCM3_xnapc_PMIL_SPRING.csv"),
  SUMMER = read.csv("Data/Seasonal/SISAL_HadCM3_xnapc_PMIL_SUMMER.csv"),
  AUTUMN = read.csv("Data/Seasonal/SISAL_HadCM3_xnapc_PMIL_AUTUMN.csv")
)

load("Data/Timeseries.RData")
