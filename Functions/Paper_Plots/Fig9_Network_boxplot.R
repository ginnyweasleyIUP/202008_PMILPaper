#################################################
## Paper Network 3 ##############################
#################################################
source("Functions/network_mask.R")
source("Functions/network_mask_ensemble.R")

#################################################

No.digits = 2
line_names = 1
col_chrono = "#C47900"
col_ensemble = "#006600"
cex_text = 1.2

#col_closest <- c("#6A51A3", "#9E9AC8", "#CBC9E2", "#F2E6F7")
#col_offset <- c("#CB181D", "#FB6A4A", "#FCAE91", "#FEE5D9")
for(plot in 1:1){
  pdf(file = paste0("Paper_Plots/Fig9_Network_boxplot.pdf"), height= PLOTTING_VARIABLES$HEIGHT, width = PLOTTING_VARIABLES$WIDTH)
  par(mar = c(5,11,1,2))
  plot(c(-1,1), c(1,28-6), type = "n", axes = FALSE, xlab = "", ylab = "" , cex = cex_text)
  abline(v=0)
  ## SITES
  corr <- list(full = list(), gauss = list(), ensemble = list())
  list <- DATA_past1000$CAVES$entity_info %>% select(site_id, entity_id) %>%
    filter(entity_id %in% DATA_past1000$CAVES$entity_info$entity_id[mask_spec]) %>% group_by(site_id) %>% count() %>% filter(n>1)
  entities <- DATA_past1000$CAVES$entity_info %>% filter(site_id %in% list$site_id) %>% 
    filter(entity_id %in% DATA_past1000$CAVES$entity_info$entity_id[mask_spec]) %>% select(entity_id)
  FULL = NETWORK$record
  FULL[NETWORK$record_p>0.1] = NA
  FULL = FULL[network_mask(list,"site")]
  corr$full = as.numeric(FULL)
  GAUSS = NETWORK$record_gauss
  GAUSS[NETWORK$record_gauss_p>0.1] = NA
  GAUSS = GAUSS[network_mask(list,"site")]
  
  corr$gauss = as.numeric(GAUSS)
  corr$ensemble = c_ensemble[network_mask_ensemble(list, "site")]
  
  boxplot(as.numeric(corr$full), add = TRUE, at = 28-6, boxwex = 1, names = "n", horizontal = T, outline = F) 
  boxplot(as.numeric(corr$gauss), add = TRUE, at = 27.5-6, boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F, outline = F)
  boxplot(as.numeric(corr$ensemble), add = TRUE, at = 27-6, boxwex = 1, names = "n", horizontal = T, col = col_ensemble, axes = F, outline = F)
  
 
  #mtext(side=2,"GMST",                cex = unitscex,    line = unitslinno, las = 1, col = "black", at = 1)
  mtext(side = 2, paste0("Site [",dim(entities)[1],"/",dim(list)[1],"]"), cex = cex_text, line = line_names, las = 1, col = "black", at = 28-6)
  abline(h=26.5-6, lty = 3, lwd = 2)
  ##GRIGBOX
  
  corr <- list(full = list(), gauss = list(), ensemble = list())
  list <- as.tibble(DATA_past1000$CAVES$gridbox_list) %>% group_by(gridbox_id) %>% count() %>% filter(n>1)
  
  entities <- as.data.frame(DATA_past1000$CAVES$gridbox_list) %>% filter(gridbox_id %in% list$gridbox_id) %>% select(entity_id)
  FULL = NETWORK$record
  FULL[NETWORK$record_p>0.1] = NA
  FULL = FULL[network_mask(list,"gridbox")]
  corr$full = as.numeric(FULL)
  GAUSS = NETWORK$record_gauss
  GAUSS[NETWORK$record_gauss_p>0.1] = NA
  GAUSS = GAUSS[network_mask(list,"gridbox")]
  corr$gauss = as.numeric(GAUSS)
  corr$gauss = as.numeric(GAUSS)
  corr$ensemble = c_ensemble[network_mask_ensemble(list, "gridbox")]
  
  boxplot(as.numeric(corr$full), add = TRUE, at = 26-6, boxwex = 1, names = "n", horizontal = T, outline = F) 
  boxplot(as.numeric(corr$gauss), add = TRUE, at = 25.5-6, boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F, outline = F)
  boxplot(as.numeric(corr$ensemble), add = TRUE, at = 25-6, boxwex = 1, names = "n", horizontal = T, col = col_ensemble, axes = F, outline = F)
  
  #mtext(side=2,"GMST",                cex = unitscex,    line = unitslinno, las = 1, col = "black", at = 1)
  mtext(side = 2, paste0("Gridbox [",dim(entities)[1],"/",dim(list)[1],"]"), cex = cex_text, line = line_names, las = 1, col = "black", at = 26-6)
  abline(h=24.5-6, lty = 3, lwd = 2)
  
  
  position <- list(cluster6 = c(24, 23.5, 23, 22.5,22)-6,
                   cluster2 = c(21, 20.5, 20, 19.5,19)-6,
                   cluster3 = c(18, 17.5, 17, 16.5,16)-6,
                   #cluster7 = c(15, 14.5, 14, 13.5,13),
                   cluster1 = c(12, 11.5, 11, 10.5,10)-3,
                   cluster5 = c( 9,  8.5,  8,  7.5,7)-3)#,
                   #cluster9 = c( 6, 5.5, 5, 4.5, 4))
  text <- list(cluster1 = "ICA", cluster2 = "S-America", cluster3 = "Europe", cluster4 = "Africa", 
               cluster5 = "CEA", cluster6 = "N-America", cluster8 = "NZ", cluster9 = "SE Asia")
  cluster_list <- as.tibble(DATA_past1000$CAVES$cluster_list) %>% group_by(cluster_id) %>% count() %>% filter(n>1)
  cluster_number <- c(6,2,3,4,7,1,5,9,8)
  
  ## CLUSTER
  #for(cluster in c(1,2,3,5,6,7)){
  for(cluster in c(6,2,3,1,5)){
    corr = list(full = list(), gauss = list(), ensemble = list(), full_sim = list(), gauss_sim = list())
    list = as.tibble(DATA_past1000$CAVES$cluster_list) %>% filter(cluster_id == cluster)
    FULL = NETWORK$record
    FULL[NETWORK$record_p>0.1] = NA
    FULL = FULL[network_mask(list,"cluster")]
    corr$full = as.numeric(FULL)
    GAUSS = NETWORK$record_gauss
    GAUSS[NETWORK$record_gauss_p>0.1] = NA
    GAUSS = GAUSS[network_mask(list,"cluster")]
    corr$gauss = as.numeric(GAUSS)
    corr$ensemble = c_ensemble[network_mask_ensemble(list, "cluster")]
    
    FULL = NETWORK$sim_full_a
    FULL[NETWORK$sim_full_a_p>0.1] = NA
    FULL = FULL[network_mask(list,"cluster")]
    corr$full_sim = as.numeric(FULL)
    GAUSS = NETWORK$sim_full_a_gauss
    GAUSS[NETWORK$sim_full_a_gauss_p>0.1] = NA
    GAUSS = GAUSS[network_mask(list,"cluster")]
    corr$gauss_sim = as.numeric(GAUSS)
    
    boxplot(corr$full, add = TRUE, at = position[[paste0("cluster",cluster)]][1] ,boxwex = 1, names = "n", horizontal = T, axes = F, outline = F) 
    boxplot(corr$gauss, add = TRUE, at = position[[paste0("cluster",cluster)]][2] , boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F, outline = F)
    boxplot(corr$ensemble, add = TRUE, at = position[[paste0("cluster",cluster)]][3] , boxwex = 1, names = "n", horizontal = T, col = col_ensemble, axes = F, outline = F)
    
    boxplot(corr$full_sim, add = TRUE, at = position[[paste0("cluster",cluster)]][4]  ,boxwex = 1, names = "n", horizontal = T, col = "dodgerblue3", axes = F, outline = F) 
    boxplot(corr$gauss_sim, add = TRUE, at = position[[paste0("cluster",cluster)]][5] , boxwex = 1, names = "n", horizontal = T, col = adjustcolor("dodgerblue3", alpha.f = 0.5), axes = F, outline = F)
    
    mtext(side = 2, paste0("c",cluster_number[cluster],"/",text[[cluster]]," [",cluster_list$n[cluster],"]"), cex = cex_text, line = line_names, las = 1, col = "black", at = position[[paste0("cluster",cluster)]][1])
    
  }
  
  abline(h=3.5, lty =3, lwd = 2)
  
  ##Global
  
  corr$full = as.numeric(NETWORK$record[NETWORK$record_p<0.1])
  corr$gauss = as.numeric(NETWORK$record_gauss[NETWORK$record_gauss_p<0.1])
  corr$full_sim = as.numeric(NETWORK$sim_full_a[NETWORK$sim_full_a_p<0.1])
  corr$gauss_sim = as.numeric(NETWORK$sim_full_a_gauss[NETWORK$sim_full_a_gauss_p<0.1])
  corr$ensembls = as.numeric(c_ensemble)
  
  boxplot(corr$full, add = TRUE, at = 3 ,boxwex = 1, names = "n", horizontal = T, axes = F) 
  boxplot(corr$gauss, add = TRUE, at = 2.5 , boxwex = 1, names = "n", horizontal = T, col = "grey", axes = F)
  boxplot(corr$ensemble, add = TRUE, at = 2 , boxwex = 1, names = "n", horizontal = T, col = col_ensemble, axes = F)
  boxplot(corr$full_sim, add = TRUE, at = 1.5,boxwex = 1, names = "n", horizontal = T, col = "dodgerblue3", axes = F, outline = F) 
  boxplot(corr$gauss_sim, add = TRUE, at = 1 , boxwex = 1, names = "n", horizontal = T, col = adjustcolor("dodgerblue3", alpha.f = 0.5), axes = F, outline = F)
  
  
  mtext(side = 2, paste0("Global [",dim(NETWORK$record)[[1]],"]"), cex = cex_text, line = line_names, las = 1, col = "black", at = 3)
  mtext(side = 1, text = "Cross-correlation between entity-pairs", cex = cex_text, col = "black", line = 2)
  
  
  legend(-1.075,-2.2,xpd = T,inset=-0.2, bty='n', x.intersp=0.5,text.width=c(1,0.3,0.5),
         c("records", "records 100yr t.sc", "records tuned"), fill = c("white", "grey", col_ensemble), horiz = TRUE, cex = cex_text)
  legend(-1.075, -3.2,xpd = T,inset=-0.2, bty='n', x.intersp=0.5,text.width=c(1,0.3),
         c("simulation", "simulation 100yr t.sc"), fill = c("dodgerblue3", adjustcolor("dodgerblue3", alpha.f = 0.5)), horiz=TRUE, cex=cex_text)
  
  
  dev.off()
  
}


rm(corr, cluster_list, position,text, cluster, cluster_number, col_chrono, col_ensemble, line_names, cex_text, plot, FULL, GAUSS, entities, No.digits, list)
