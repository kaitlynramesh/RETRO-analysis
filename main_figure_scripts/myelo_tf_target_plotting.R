# Visualize TF/target gene pairs + simulation
# Figure 4F-I

library(ggpubr) # plotting
source("~/KaitlynRRStudio/RETRO-analysis/main_figure_scripts/tf_modeling_functions.R")

# relevant functions
fitting_gene_pair <- function(gene_pair, curves, init_param, lag=0, upper_bounds=c(100, 100, 2, 5, 10)) {
  
  # Obtain TF and target trajectories from PseudotimeDE
  curves_tf = curves[["TF"]]
  curves_tar = curves[["target"]]
  
  # Obtain TF and target trajectories from PseudotimeDE
  tf_fitting <- curves_tf[["data"]][curves_tf[["data"]]$gene==gene_pair[1],c(2,4)]
  target_fitting <- curves_tar[["data"]][curves_tar[["data"]]$gene==gene_pair[2],c(2,4)]
  
  t <- tf_fitting$pseudotime # pseudotime values
  
  # Generate TF v. target data frame
  if(lag!=0) {
    x = unlist(tf_fitting[order(t),2])
    y = unlist(target_fitting[order(t),2])
    
    t = sort(t) # make sure pseudotime is ordered
    shift = which.min(abs(t-lag)) # how many pseudotime-points over
    
    x_lag = x[1:(length(x) - shift)] # LAGGED TF expression
    y_lag = y[(shift + 1):length(y)] # LAGGED target expression
    
    t_early = t[1:(length(t) - shift)] # shifted time trajectories
    t_late = t[(shift + 1):length(t)] 
    
    tf_vs_target = data.frame(cbind(t_late, x_lag, y_lag))
    
  } else {
    tf_vs_target <- data.frame(cbind(tf_fitting[order(t), ], 
                                     target_fitting$fitted[order(target_fitting$pseudotime)]))
    t_early = sort(t)
    t_late = sort(t)
  }
  colnames(tf_vs_target) <- c('Pseudotime', 'V2', 'V3')
  
  # Optimization
  lower_range = c(0,0,k=1e-5,0,0)
  upper_range = upper_bounds
  results_simu <- model_sim(init_param, tf_vs_target, target_model, 
                            lower_range, upper_range)
  param_est = results_simu[[3]]
  
  # Simulation parameters
  t = tf_vs_target$Pseudotime
  X = tf_vs_target[,2] # TF
  Xn = tf_vs_target[1,2:3] # (TF, target)
  t.total = t
  dt = diff(t)
  sim_target = RK4_generic(derivs = target_model, Xn = Xn, X = X,
                           t.total = t.total, dt = dt, param = param_est)
  
  # Plot trajectories
  plot(t_early, tf_vs_target[,2], ylim=c(0,6), type='l', xlab='Time', ylab='Counts')
  points(t_late, tf_vs_target[,3], type='l', col='blue')
  points(t_late, sim_target, col='red', type='l')
  legend('topright', c(gene_pair[1], gene_pair[2], paste0(gene_pair[2], ' (sim)')), cex=.45, 
         col=c('black', 'blue', 'red'), lty=1)
  
  
  simulation_res = data.frame(cbind(pseudotime_tf = t_early, 
                                    pseudotime_target = t_late,
                                    tf = X,
                                    target = tf_vs_target[,3],
                                    simulation = sim_target))
  return(simulation_res)
}


#### (1) Megakaryocyte Lineage #### 
mega_gene_vec = c(TF = "GTF2I", target = "HSPA5")
gene_class = c("TF", "target")
res_1 = tf_target_comp_meg[["results"]] # pseudotimeDE param fitting
index_1 = match(mega_gene_vec, res_1$gene)  # indexing

curves_1 = vector(mode="list", length=length(mega_gene_vec))
for(i in seq(index_1)) {
  
  ind = index_1[i]
  gene_curve <- PseudotimeDE::plotCurve(gene.vec = res_1$gene[ind],
                                        ori.tbl = retro_tbl_meg[[2]],
                                        assay = 'logcounts',
                                        mat = retro_tbl_meg[[1]],
                                        model.fit = res_1$gam.fit[ind])
  
  gene_header_col = ifelse(gene_class[i] == "TF", "#F89C86", "#95A3F3")
  gene_curve[["theme"]][["strip.background"]][["fill"]] = gene_header_col
  
  curves_1[[i]] = gene_curve
}

curves_1
names(curves_1) = gene_class

#### (2) Macrophage Lineage #### 
mac_gene_vec = c(TF = "CEBPB", target = "ITGB2")
gene_class = c("TF", "target")
res_2 = tf_target_comp_mac[["results"]] # pseudotimeDE param fitting
index_2 = match(mac_gene_vec, res_2$gene)  # indexing

curves_2 = vector(mode="list", length=length(mac_gene_vec))
for(i in seq(index_2)) {
  
  ind = index_2[i]
  gene_curve <- PseudotimeDE::plotCurve(gene.vec = res_2$gene[ind],
                                        ori.tbl = retro_tbl_mac[[2]],
                                        assay = 'logcounts',
                                        mat = retro_tbl_mac[[1]],
                                        model.fit = res_2$gam.fit[ind])
  
  gene_header_col = ifelse(gene_class[i] == "TF", "#F89C86", "#95A3F3")
  gene_curve[["theme"]][["strip.background"]][["fill"]] = gene_header_col
  
  curves_2[[i]] = gene_curve
}

curves_2
names(curves_2) = gene_class

#### (3) Monocyte Lineage ####
mon_gene_vec = c(TF = "SPI1", target="S100A9")
gene_class = c("TF", "target")
res_3 = tf_target_comp_mon[["results"]] # pseudotimeDE param fitting
index_3 = match(mon_gene_vec, res_3$gene)  # indexing

curves_3 = vector(mode="list", length=length(mon_gene_vec))
for(i in seq(index_3)) {
  
  ind = index_3[i]
  gene_curve <- PseudotimeDE::plotCurve(gene.vec = res_3$gene[ind],
                                        ori.tbl = retro_tbl_mon[[2]],
                                        assay = 'logcounts',
                                        mat = retro_tbl_mon[[1]],
                                        model.fit = res_3$gam.fit[ind])
  
  gene_header_col = ifelse(gene_class[i] == "TF", "#F89C86", "#95A3F3")
  gene_curve[["theme"]][["strip.background"]][["fill"]] = gene_header_col
  
  curves_3[[i]] = gene_curve
}

curves_3
names(curves_3) = gene_class


#### (4) Progenitor Cell State ####
mast_gene_vec = c(TF = "POU5F1", target = "TDGF1")
gene_class = c("TF", "target")
res_4 = tf_target_comp_mast[["results"]] # pseudotimeDE param fitting
index_4 = match(mast_gene_vec, res_4$gene)  # indexing

curves_4 = vector(mode="list", length=length(mast_gene_vec))
for(i in seq(index_4)) {
  
  ind = index_4[i]
  gene_curve <- PseudotimeDE::plotCurve(gene.vec = res_4$gene[ind],
                                        ori.tbl = retro_tbl_mast[[2]],
                                        assay = 'logcounts',
                                        mat = retro_tbl_mast[[1]],
                                        model.fit = res_4$gam.fit[ind])
  
  gene_header_col = ifelse(gene_class[i] == "TF", "#F89C86", "#95A3F3")
  gene_curve[["theme"]][["strip.background"]][["fill"]] = gene_header_col
  
  curves_4[[i]] = gene_curve
}

curves_4
names(curves_4) = gene_class

(curves_4[[1]] + curves_4[[2]])  # progenitor
(curves_1[[1]] + curves_1[[2]]) # mega

(curves_3[[1]] + curves_3[[2]])  # mono
(curves_2[[1]] + curves_2[[2]]) # mac


#### TF/Target Fitting + Analysis ####

# Macrophage fitting
mac_gene_vec = c("CEBPB", "ITGB2")
init_param_1 = c(g0=0.00404862, lam=65.4432, k=0.374534, n=10.0000, th=0.323116)
sim_res_1 = fitting_gene_pair(mac_gene_vec, curves_2, init_param_1, lag=0,
                              upper_bounds = c(100,100,2,10,1))
g1 = ggplot(data=sim_res_1) +  
  geom_line(aes(x=pseudotime_tf, y=tf, color = mac_gene_vec[1]), lwd=1) +
  geom_line(aes(x=pseudotime_target, y=target, color = mac_gene_vec[2]), lwd=1) + 
  geom_line(aes(x=pseudotime_target, y=simulation, color = paste0(mac_gene_vec[2], " (sim)")), lwd=1) +
  scale_color_manual(
    name = "gene",
    values = setNames(
      c("#E64B35", "#8C6BB1", "#3C4EC2"),
      c(mac_gene_vec[1], mac_gene_vec[2], paste0(mac_gene_vec[2], " (sim)"))
    )
  ) + theme_bw() + 
  geom_vline(xintercept=31.1279, color = "red", linetype = "dashed") + 
  ylab("log(counts+1)") + 
  xlab("pseudotime") + 
  xlim(min(pseudotime), max(pseudotime)) + 
  guides(colour="none")

# p_est = 0.00404864  65.4432 0.374534  10.0000 0.323116


##

# Monocyte fitting
init_param_2 = c(g0=0.004, lam=65, k=0.3, n=10, th=0.3)
sim_res_2 = fitting_gene_pair(mon_gene_vec, curves_3, init_param_2, lag=1.263158, 
                              upper_bound =c(100,100,2,10,10))
g2 = ggplot(data=sim_res_2) +  
  geom_line(aes(x=pseudotime_tf, y=tf, color = mon_gene_vec[1]), lwd=1) +
  geom_line(aes(x=pseudotime_target, y=target, color = mon_gene_vec[2]), lwd=1) + 
  geom_line(aes(x=pseudotime_target, y=simulation, color = paste0(mon_gene_vec[2], " (sim)")), lwd=1) +
  scale_color_manual(
    name = "gene",
    values = setNames(
      c("#E64B35", "#8C6BB1", "#3C4EC2"),
      c(mon_gene_vec[1], mon_gene_vec[2], paste0(mon_gene_vec[2], " (sim)"))
    )
  ) + theme_bw() + 
  geom_vline(xintercept=22.9229421, color = "red", linetype = "dashed") + 
  ylab("log(counts+1)") + 
  xlab("pseudotime") +
  xlim(min(pseudotime), max(pseudotime)) + 
  guides(colour="none")

# p_Est = 0.00811769  28.2806 0.0947161  1.31836 0.000236397


## 

# Mast cell fitting
mast_gene_vec = c("POU5F1", "TDGF1")
init_param_3 = c(g0=.04, lam=300, k=1, n=1, th=0.3)
sim_res_3 = fitting_gene_pair(mast_gene_vec, curves_4, init_param_3, lag=0) #, upper_bounds = c(100,400,1,5,10))
g3 = ggplot(data=sim_res_3) +  
  geom_line(aes(x=pseudotime_tf, y=tf, color = mast_gene_vec[1]), lwd=1) +
  geom_line(aes(x=pseudotime_target, y=target, color = mast_gene_vec[2]), lwd=1) + 
  geom_line(aes(x=pseudotime_target, y=simulation, color = paste0(mast_gene_vec[2], " (sim)")), lwd=1) +
  scale_color_manual(
    name = "gene",
    values = setNames(
      c("#E64B35", "#8C6BB1", "#3C4EC2"),
      c(mast_gene_vec[1], mast_gene_vec[2], paste0(mast_gene_vec[2], " (sim)"))
    )
  ) + theme_bw() + 
  geom_vline(xintercept=0.6713357, color = "red", linetype = "dashed") + 
  ylab("log(counts+1)") + 
  xlab("pseudotime") + 
  xlim(min(pseudotime), max(pseudotime)) +
  guides(colour="none")

# p_Est = 7.27700e-08  99.9958 0.161178  1.48568  3.89490


# Megakaryocyte fitting
gene_pair_4 = c("GTF2I", "HSPA5")
init_param_4 = c(g0=0.035, lam=100, k=0.23, n=2.15, th=2.7)
sim_res_4 = fitting_gene_pair(gene_pair_4, curves_1, init_param_4, lag=0, upper_bound =c(100,200,2,5,10))
g4 = ggplot(data=sim_res_4, aes(x=pseudotime)) +  
  geom_line(aes(x=pseudotime_tf, y=tf, color = gene_pair_4[1]), lwd=1) +
  geom_line(aes(x=pseudotime_target, y=target, color = gene_pair_4[2]), lwd=1) + 
  geom_line(aes(x=pseudotime_target, y=simulation, color = paste0(gene_pair_4[2], " (sim)")), lwd=1) +
  scale_color_manual(
    name = "gene",
    values = setNames(
      c("#E64B35", "#8C6BB1", "#3C4EC2"),
      c(gene_pair_4[1], gene_pair_4[2], paste0(gene_pair_4[2], " (sim)"))
    )
  ) +  
  theme_bw() + 
  geom_vline(xintercept=8.2663187, color = "red", linetype = "dashed") + 
  ylab("log(counts+1)") +
  xlab("pseudotime") + 
  xlim(min(pseudotime), max(pseudotime)) + 
  guides(colour="none")

# p_Est =  0.153231  112.908 0.301709  5.00000  1.94872

ggarrange(g3,g4,g2,g1, ncol=1)




