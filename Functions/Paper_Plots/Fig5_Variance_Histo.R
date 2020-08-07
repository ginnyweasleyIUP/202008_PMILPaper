
#################################################

run = "a"

# Variances for simulation and for proxy:

VARIANCE <- list(
  var_d18Oc = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_proxy_a = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_proxy_b = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_proxy_c = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_a = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_b = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_c = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_a = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_b = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_c = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_a_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_b_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_ds_c_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_a_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_b_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id)),
  var_full_c_isot = numeric(length(DATA_past1000$CAVES$entity_info$entity_id))
)

for(ii in 1:length(DATA_past1000$CAVES$entity_info$entity_id)){
  entity = DATA_past1000$CAVES$entity_info$entity_id[ii]
  data_rec = DATA_past1000$CAVES$record_res %>% filter(entity_id == entity)
  data_yearly_a = DATA_past1000$CAVES$yearly_res$a %>% filter(entity_id == entity)
  data_yearly_b = DATA_past1000$CAVES$yearly_res$b %>% filter(entity_id == entity)
  data_yearly_c = DATA_past1000$CAVES$yearly_res$c %>% filter(entity_id == entity)
  VARIANCE$var_d18Oc[ii] = var(data_rec$d18O_measurement, na.rm = T)
  VARIANCE$var_proxy_a[ii] = var(data_rec$d18O_dw_eq_a, na.rm = T)
  VARIANCE$var_proxy_b[ii] = var(data_rec$d18O_dw_eq_b, na.rm = T)
  VARIANCE$var_proxy_c[ii] = var(data_rec$d18O_dw_eq_c, na.rm = T)
  VARIANCE$var_ds_a[ii] = var(data_rec$ITPC_a, na.rm = T)
  VARIANCE$var_ds_b[ii] = var(data_rec$ITPC_b, na.rm = T)
  VARIANCE$var_ds_c[ii] = var(data_rec$ITPC_c, na.rm = T)
  VARIANCE$var_full_a[ii] = var(data_yearly_a$ITPC, na.rm = T)
  VARIANCE$var_full_b[ii] = var(data_yearly_b$ITPC, na.rm = T)
  VARIANCE$var_full_c[ii] = var(data_yearly_c$ITPC, na.rm = T)
  VARIANCE$var_ds_a_isot[ii] = var(data_rec$ISOT_a, na.rm = T)
  VARIANCE$var_ds_b_isot[ii] = var(data_rec$ISOT_b, na.rm = T)
  VARIANCE$var_ds_c_isot[ii] = var(data_rec$ISOT_c, na.rm = T)
  VARIANCE$var_full_a_isot[ii] = var(data_yearly_a$ISOT, na.rm = T)
  VARIANCE$var_full_b_isot[ii] = var(data_yearly_b$ISOT, na.rm = T)
  VARIANCE$var_full_c_isot[ii] = var(data_yearly_c$ISOT, na.rm = T)
}

VARIANCE$var_proxy <- as.numeric(c(VARIANCE$var_proxy_a, VARIANCE$var_proxy_b, VARIANCE$var_proxy_c))
VARIANCE$var_ds <- as.numeric(c(VARIANCE$var_ds_a, VARIANCE$var_ds_b, VARIANCE$var_ds_c))
VARIANCE$var_full <- as.numeric(c(VARIANCE$var_full_a, VARIANCE$var_full_b, VARIANCE$var_full_c))
VARIANCE$var_ds_isot <- as.numeric(c(VARIANCE$var_ds_a_isot, VARIANCE$var_ds_b_isot, VARIANCE$var_ds_c_isot))
VARIANCE$var_full_isot <- as.numeric(c(VARIANCE$var_full_a_isot, VARIANCE$var_full_b_isot, VARIANCE$var_full_c_isot))


## PLOTTING

cex_text = 1.7

#lines(c(median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS$VR_sim_temp)),median(log10(PLOTTING_DATA$FIG4$VAR_RATIOS$VR_sim_temp))), c(0, ))
  xlimz = c(-2.2,2.2)
  title_pos = -2.2
  pdf(file = paste0("Paper_Plots/Fig5_Variance_2_histo_.pdf"), width = 1.3*6, height = 2*PLOTTING_VARIABLES$HEIGHT/1.5)
  #png(file = paste0("Plots/Paper_Plot_4_Variance_2_histo_xnap",run,".png"), width = 30*2*6, height = 50*2*PLOTTING_VARIABLES$HEIGHT/1.5)
  par(mfrow=c(2,1),oma = c(1,3,0,0) + 0.1,mar = c(3,1,0,1) + 0.1, new = FALSE)

  hist(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full_isot")]][rep(mask_var,3)]),
       breaks = 9, border = "white", prob = TRUE,
       ylim = c(0,1), xlim = xlimz, xlab = "",xaxt = 'n',
       main = "", cex.main = cex_text, cex.axis = cex_text)
  axis(side = 1, at = c(log10(0.01), log10(0.1), 0, log10(10), log10(100)),
       labels = c(0.01, 0.1, 1, 10, 100), cex.axis = cex_text)
  lines(c(median(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full_isot")]][rep(mask_var,3)]), na.rm = T),
          median(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full_isot")]][rep(mask_var,3)]), na.rm = T)),
        c(0, max(density(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full_isot")]][rep(mask_var,3)]), na.rm = T)$y)-0.01),
        lwd = 2, col = "black", lty = 2)
  lines(density(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full_isot")]][rep(mask_var,3)]), na.rm = T),
        lwd = 2, col = "black")
  lines(c(median(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds_isot")]][rep(mask_var,3)]), na.rm = T),
          median(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds_isot")]][rep(mask_var,3)]), na.rm = T)),
        c(0, max(density(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds_isot")]][rep(mask_var,3)]), na.rm = T)$y)-0.01),
        lwd = 2, col = "#B2182B", lty = 2)
  lines(density(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds_isot")]][rep(mask_var,3)]),na.rm = T),
        lwd = 2, col = "#B2182B")
  abline(v=0, col = "grey60", lty = 3)
  #mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.5, cex = 1.5)
  mtext(text = "density",side = 2,line = 2.5, cex = cex_text)
  text(1.0, 0.8, "down-sampled", col = "#B2182B", cex = cex_text)
  text(-1.4, 0.6, "full", col = "black", cex = cex_text)
  mtext(text = expression(paste(delta^{plain(18)}, plain(O))), side = 3, line = -2, adj = 0,col = "black", cex = cex_text, at = title_pos)


  hist(log10(VARIANCE[[paste0("var_proxy_", run)]][rep(mask_var,3)]/VARIANCE$var_full_a[rep(mask_var,3)]),
       breaks = 9, border = "white", prob = TRUE,
       ylim = c(0,1), xlim = xlimz, xlab = "",xaxt = 'n',
       main = "", cex.main = cex_text, cex.axis = cex_text)
  axis(side = 1, at = c(log10(0.01), log10(0.1), 0, log10(10), log10(100)),
       labels = c(0.01, 0.1, 1, 10, 100), cex.axis = cex_text)
  lines(c(median(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full")]][rep(mask_var,3)]), na.rm = T),
          median(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full")]][rep(mask_var,3)]), na.rm = T)),
        c(0, max(density(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full")]][rep(mask_var,3)]), na.rm = T)$y)-0.01),
        lwd = 2, col = "black", lty = 2)
  lines(density(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_full")]][rep(mask_var,3)]), na.rm = T),
        lwd = 2, col = "black")
  lines(c(median(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds")]][rep(mask_var,3)]), na.rm = T),
          median(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds")]][rep(mask_var,3)]), na.rm = T)),
        c(0, max(density(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds")]][rep(mask_var,3)]), na.rm = T)$y)-0.05),
        lwd = 2, col = "#B2182B", lty = 2)
  lines(density(log10(VARIANCE[[paste0("var_proxy")]][rep(mask_var,3)]/VARIANCE[[paste0("var_ds")]][rep(mask_var,3)]), na.rm = T),
        lwd = 2, col = "#B2182B")
  abline(v=0, col = "grey60", lty = 3)
  mtext(text = "density",side = 2,line = 2.5, cex = cex_text)
  mtext(text = expression(paste(delta^{plain(18)}, plain(O)[plain(pw)])), side = 3, line = -2, adj = 0,col = "black", cex = cex_text, at = title_pos)
  mtext(text = TeX("$$Var_{Rec}/Var_{Sim}$$"),side = 1,line = 2.7, cex = cex_text)

  dev.off()



remove(ii, entity, run, xlimz, VARIANCE, data_rec, data_yearly_a, data_yearly_b, data_yearly_c, title_pos, cex_text)


