#################################################
## DISCUSSION NETWORK TUNING ####################
#################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(zoo)
library(PaleoSpec)
library(nest)
library(latex2exp)

NETWORK_TUNING <- list()

ensemble_entities <- c(14,21,33,48,51,76,85,93,94,95,97,111,113,117,118,123,128,137,147,165,172,178,179,185,187,
                       201,202,209,212,222,226,238,240,242,253,278,286,289,294,296,305,326,329,330,335,351,361,378,390,395,399,
                       420,422,430,435,437,442,443,447,448,461,464,496,498,506,514,522,528,538,539,541,544,546,547,548,560,564,573,577,588,589,595,598,
                       613,620,621,623,672,673)

C_ensemble_new <- cbind(c_ensemble[,1:6],c_ensemble[,8:11], c_ensemble[,13:89])
C_ensemble_new <- rbind(C_ensemble_new[1:6,], C_ensemble_new[8:11,], C_ensemble_new[13:89,])
C_fill_gaps <- NETWORK$record
C_fill_gaps[NETWORK$record_p>0.1] = NA
C_ensemble_new[is.na(C_ensemble_new)] = C_fill_gaps[is.na(C_ensemble_new)]

#################################################

# SITE
NETWORK_TUNING$SITES <- list()
list <- DATA_past1000$CAVES$entity_info %>% select(site_id, entity_id) %>%
  filter(entity_id %in% DATA_past1000$CAVES$entity_info$entity_id[mask_spec]) %>% group_by(site_id) %>% count() %>% filter(n>1)
entities <- DATA_past1000$CAVES$entity_info %>% filter(site_id %in% list$site_id) %>% 
  filter(entity_id %in% DATA_past1000$CAVES$entity_info$entity_id[mask_spec]) %>% select(entity_id)
FULL = NETWORK$record
FULL[NETWORK$record_p>0.1] = NA
FULL[!network_mask(list,"site")] = NA
NETWORK_TUNING$SITES$C_rec = FULL
GAUSS = NETWORK$record_gauss
GAUSS[NETWORK$record_gauss_p>0.1] = NA
GAUSS[!network_mask(list,"site")] = NA
NETWORK_TUNING$SITES$C_rec_gauss = GAUSS
ENSEMBLE = C_ensemble_new
ENSEMBLE[!network_mask(list, "site")] = NA
NETWORK_TUNING$SITES$C_rec_ensemble = ENSEMBLE

# GRIDBOX 

NETWORK_TUNING$GRIDBOX <- list()
list <- as.tibble(DATA_past1000$CAVES$gridbox_list) %>% group_by(gridbox_id) %>% count() %>% filter(n>1)
FULL = NETWORK$record
FULL[NETWORK$record_p>0.1] = NA
FULL[!network_mask(list,"gridbox")] = NA
NETWORK_TUNING$GRIDBOX$C_rec = FULL
GAUSS = NETWORK$record_gauss
GAUSS[NETWORK$record_gauss_p>0.1] = NA
GAUSS[!network_mask(list,"gridbox")] = NA
NETWORK_TUNING$GRIDBOX$C_rec_gauss = GAUSS
ENSEMBLE = C_ensemble_new
ENSEMBLE[!network_mask(list, "gridbox")] = NA
NETWORK_TUNING$GRIDBOX$C_rec_ensemble = ENSEMBLE

# CLUSTER

NETWORK_TUNING$CLUSTER <- list()

for(cluster in 1:9){
  
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]] <- list()
  list = as.tibble(DATA_past1000$CAVES$cluster_list) %>% filter(cluster_id == cluster)
  
  FULL = NETWORK$record
  FULL[NETWORK$record_p>0.1] = NA
  FULL[!network_mask(list, "cluster")] = NA
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]]$C_rec = FULL
  GAUSS = NETWORK$record_gauss
  GAUSS[NETWORK$record_gauss_p>0.1] = NA
  GAUSS[!network_mask(list,"cluster")] = NA
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]]$C_rec_gauss = FULL
  ENSEMBLE = C_ensemble_new
  ENSEMBLE[!network_mask(list, "cluster")] = NA
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]]$C_rec_ensemble = ENSEMBLE
  
  ## mask 50% closest within cluster
  dist_mask_small <- matrix(logical(length = sum(mask_spec)^2), ncol = sum(mask_spec))
  dist_mask_small[lower.tri(dist_mask_small)] = fossil::earth.dist(cbind(DATA_past1000$CAVES$entity_info$latitude[mask_spec],DATA_past1000$CAVES$entity_info$longitude[mask_spec]),dist=TRUE)
  dist_mask_small[dist_mask_small == 0] = NA
  dist_mask_small[!network_mask(list, "cluster")] = NA
  dist_mask_small[dist_mask_small>median(dist_mask_small, na.rm = T)] = NA
  
  mask_small = matrix(logical(length = sum(mask_spec)^2), ncol = sum(mask_spec))
  mask_small[!is.na(dist_mask_small)] = T
  
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest <- list()
  FULL = NETWORK$record
  FULL[NETWORK$record_p>0.1] = NA
  FULL[!mask_small] = NA
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]]$closest$C_rec = FULL
  GAUSS = NETWORK$record_gauss
  GAUSS[NETWORK$record_gauss_p>0.1] = NA
  GAUSS[!mask_small] = NA
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]]$closest$C_rec_gauss = FULL
  ENSEMBLE = C_ensemble_new
  ENSEMBLE[!mask_small] = NA
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]]$closest$C_rec_ensemble = ENSEMBLE
  
  ## mask 50% smallest offset
  offset <- numeric(length(list$entity_id))
  for(ii in 1:length(list$entity_id)){
    entity = list$entity_id[ii]
    data_rec <- DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
    offset[ii] = mean(data_rec$ITPC_a, na.rm = T) - mean(data_rec$d18O_dw_eq_a, na.rm = T)
  }
  
  list_new = list
  list_new$entity_id[abs(offset)>median(abs(offset), na.rm = T)] = NA
  list_new = na.omit(list_new)
  
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset <- list()
  FULL = NETWORK$record
  FULL[NETWORK$record_p>0.1] = NA
  FULL[!network_mask(list_new, "cluster")] = NA
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]]$smallest_offset$C_rec = FULL
  GAUSS = NETWORK$record_gauss
  GAUSS[NETWORK$record_gauss_p>0.1] = NA
  GAUSS[!network_mask(list_new,"cluster")] = NA
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]]$smallest_offset$C_rec_gauss = FULL
  ENSEMBLE = C_ensemble_new
  ENSEMBLE[!network_mask(list_new, "cluster")] = NA
  NETWORK_TUNING$CLUSTER[[paste0("CLUSTER", cluster)]]$smallest_offset$C_rec_ensemble = ENSEMBLE
  
}

#################################################
## EXTENDED BOXPLOT #############################
#################################################

#sites = 4
#gridbox = 4
#cluster = 7*6 = 42
#global = 6
 
#total = 56+3 = 59

line_names = 1
col_ensemble = "#004F00"
cex_text = 1

col_closest <- c("#6A51A3", "#9E9AC8", "#CBC9E2", "#F2E6F7")
col_offset <- c("#CB181D", "#FB6A4A", "#FCAE91", "#FEE5D9")
for(plot in 1:1){
  pdf(file = paste0("Sup_Plots/SF_Network_boxplot_extended.pdf"), height= PLOTTING_VARIABLES$HEIGHT, width = 1.5*PLOTTING_VARIABLES$WIDTH)
  par(mar = c(5,11,1,2))
  plot(c(-1,1), c(1,24), type = "n", axes = FALSE, xlab = "", ylab = "" )
  abline(v=0)
  position <- list(cluster6 = c(19, 18.5, 18, 17.5, 17, 16.5, 16)+5,
                   cluster2 = c(15, 14.5, 14, 13.5, 13, 12.5, 12)+4,
                   cluster3 = c(11, 10.5, 10,  9.5,  9,  8.5, 8)+3,
                   cluster1 = c( 7,  6.5,  6,  5.5,  5,  4.5, 4)+2,
                   cluster5 = c( 4,  3.5,  3,  2.5,  2,  1.5, 1))
  text <- list(cluster1 = "India", cluster2 = "SouthAm", cluster3 = "Europe", cluster4 = "Africa", 
               cluster5 = "China", cluster6 = "NorthAm", cluster7 = "Arabia", cluster8 = "NZ", cluster9 = "SE Asia")
  cluster_list <- as.tibble(DATA_past1000$CAVES$cluster_list) %>% group_by(cluster_id) %>% count() %>% filter(n>1)
  cluster_number <- c(6,2,3,4,7,1,5,9,8)
  
  ## CLUSTER
  
  #for(cluster in c(1,2,3,6,7)){
  for(cluster in  c(6,2,3,1,5)){
    
    if(cluster == 6){
      boxplot(as.numeric(NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$C_rec), add = TRUE, at = position[[paste0("cluster",cluster)]][1] ,boxwex = 1, names = "n", horizontal = T, outline = F) 
    }else{
      boxplot(as.numeric(NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$C_rec), add = TRUE, at = position[[paste0("cluster",cluster)]][1] ,boxwex = 1, names = "n", horizontal = T, axes = F, outline = F) 
    }
    
    boxplot(as.numeric(NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest$C_rec), add = TRUE, at = position[[paste0("cluster",cluster)]][2] , boxwex = 1, names = "n", horizontal = T, col = col_closest[1], axes = F, outline = F)
    boxplot(as.numeric(NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest$C_rec_gauss), add = TRUE, at = position[[paste0("cluster",cluster)]][3] , boxwex = 1, names = "n", horizontal = T, col = col_closest[2], axes = F, outline = F)
    boxplot(as.numeric(NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$closest$C_rec_ensemble), add = TRUE, at = position[[paste0("cluster",cluster)]][4] , boxwex = 1, names = "n", horizontal = T, col = col_closest[4], axes = F, outline = F)
    
    boxplot(as.numeric(NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset$C_rec), add = TRUE, at = position[[paste0("cluster",cluster)]][5] , boxwex = 1, names = "n", horizontal = T, col = col_offset[1], axes = F, outline = F)
    boxplot(as.numeric(NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset$C_rec_gauss), add = TRUE, at = position[[paste0("cluster",cluster)]][6] , boxwex = 1, names = "n", horizontal = T, col = col_offset[2], axes = F, outline = F)
    boxplot(as.numeric(NETWORK_TUNING$CLUSTER[[paste0("CLUSTER",cluster)]]$smallest_offset$C_rec_ensemble), add = TRUE, at = position[[paste0("cluster",cluster)]][7] , boxwex = 1, names = "n", horizontal = T, col = col_offset[4], axes = F, outline = F)
    
    mtext(side = 2, paste0("c",cluster_number[cluster],"/",text[[cluster]]," [",cluster_list$n[cluster],"]"), cex = 1, line = line_names, las = 1, col = "black", at = position[[paste0("cluster",cluster)]][1]-1)
    
  }
  
  mtext(side = 1, text = "Cross-correlation between entity-pairs", cex = cex_text, col = "black", line = 1.7)
  
  legend(-1.075,-2.3,xpd = T,inset=-0.2, bty='n', x.intersp=0.5,text.width=c(1,0.25,0.4, 0.37),
         c("record", "closest record","cl.record 100yr t.sc","closest record ensemble"), fill=c("white", col_closest[1],col_closest[2],col_closest[4]), horiz=TRUE, cex=1)
  legend(-1.075,-3.3,xpd = T,inset=-0.2, bty='n', x.intersp=0.5,text.width=c(1,0.25,0.4, 0.37),
         c("offset record","offset record 100yr t.sc","offset record ensemble"), fill=col_offset, horiz=TRUE, cex=1)
  # 
  # legend(-1.075,-2.2,xpd = T,inset=-0.2, bty='n', x.intersp=0.5,text.width=c(1,0.3,0.5),
  #        c("records", "records 100yr t.sc", "records tuned"), fill = c("white", "grey", col_ensemble), horiz = TRUE, cex = cex_text)
  # legend(-1.075, -3.2,xpd = T,inset=-0.2, bty='n', x.intersp=0.5,text.width=c(1,0.3),
  #        c("simulation", "simulation 100yr t.sc"), fill = c("dodgerblue3", adjustcolor("dodgerblue3", alpha.f = 0.5)), horiz=TRUE, cex=cex_text)

  
  dev.off()
  
}

rm(cex_text, cluster, col_closest, position, mask_small, list, list_new, FULL, GAUSS, ENSEMBLE, dist_mask_small, data_rec, C_fill_gaps, C_ensemble_new)
rm(ensemble_entities, entity, ii, line_names, offset, plot, text)
rm(col_ensemble, col_offset, cluster_number, entities, cluster_list)
