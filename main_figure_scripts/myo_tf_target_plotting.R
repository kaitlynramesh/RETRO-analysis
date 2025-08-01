

#### TF fitting and simulation functions ####
source("~/KaitlynRRStudio/PseudotimeProject/TF_Modeling_Functions.R")

# Simulates target gene expression given lagged TF trajectory
fitting_gene_pair <- function(gene_pair, curves, init_param, lag=0, upper_bounds=c(100, 100, 2, 5, 10)) {
  
  # Obtain TF and target trajectories from PseudotimeDE
  tf_fitting <- curves[["data"]][curves[["data"]]$gene==gene_pair[1],c(2,4)]
  target_fitting <- curves[["data"]][curves[["data"]]$gene==gene_pair[2],c(2,4)]
  
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


#### Load datasets and pseudotime object ####
myo_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/real/scd_myo.rds"
scd_myo = readRDS(gzcon(url(myo_url)))
time = scd_myo@phenoData@data[["time"]]

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
retro_pt_obj = readRDS(file=paste0(dir, "RETRO_PT_myo.rds"))
pseudotime = retro_pt_obj@pseudotime

# pre-processed MI data
load("~/KaitlynRRStudio/RETRO-analysis/tf_target_data/myo_retro_tbl.rda")


#### Run PseudotimeDE to get trajectories 
gene_vec = c(c("Fos", "Egr1", "Jun", "Thra", "Nfia", "Ar", "Tgfb1", "Ctnnb1", "Tcf4", "Jund", "Zeb1"),
             c("Mt1", "Cd9", "Cxcl1", "Col6a1", "Rbp1", "Lpl", "Col1a2", "Dpep1", "Emp1", "Wisp2", "Acta2", "Mmp2", "Tagln2"))
knots = c(5, 10, 15) 

res <- PseudotimeDE::runPseudotimeDE(gene.vec = gene_vec,
                                     ori.tbl = retro_tbl[[2]],
                                     sub.tbl = retro_tbl[[3]], 
                                     mat = retro_tbl[[1]], ## You can also use a matrix or SeuratObj as input
                                     mc.cores = 4,
                                     assay.use = "logcounts",
                                     model = "nb", # bc log-transformed
                                     knots=knots,
                                     k=length(knots))

curves <- PseudotimeDE::plotCurve(gene.vec = res$gene,
                                  ori.tbl = retro_tbl[[2]],
                                  assay = 'logcounts',
                                  mat = retro_tbl[[1]],
                                  model.fit = res$gam.fit,
                                  ncol=5)
curves



#### Simulation process ####

gene_pair_1 = c("Jun", "Cxcl1") 
init_param_1 = c(g0=0.3, lam=90, k=0.8, th=2.7, n=8)
sim_res_1 = fitting_gene_pair(gene_pair_1, curves, init_param_1, lag=0,
                              upper_bounds = c(100,100,1,5,5))
# plotting
g1 = ggplot(data=sim_res_1) +  
  geom_line(aes(x=pseudotime_tf, y=tf, color = gene_pair_1[1]), lwd=0.75) +
  geom_line(aes(x=pseudotime_target, y=target, color = gene_pair_1[2]), lwd=0.75) + 
  geom_line(aes(x=pseudotime_target, y=simulation, color = paste0(gene_pair_1[2], " (sim)")), lwd=0.75) +
  scale_color_manual(
    name = "gene",
    values = setNames(
      c("#E64B35", "#95A3F3", "#0D1B6A"),
      c(gene_pair_1[1], gene_pair_1[2], paste0(gene_pair_1[2], " (sim)"))
    )
  ) + theme_bw() + 
  ylab("log10(counts+1)") + 
  xlab("pseudotime") + 
  guides(colour="none")


gene_pair_2 = c("Nfia", "Rbp1")
init_param_2 = c(g0=0.304515, lam=201, k=4.72708, n=9.54024, th=1.30692)
sim_res_2 = fitting_gene_pair(gene_pair_2, curves, init_param_2, lag=0,
                              upper_bounds = c(500,500,5,10,5))

g3 = ggplot(data=sim_res_2) +  
  geom_line(aes(x=pseudotime_tf, y=tf, color = gene_pair_2[1]), lwd=0.75) +
  geom_line(aes(x=pseudotime_target, y=target, color = gene_pair_2[2]), lwd=0.75) + 
  geom_line(aes(x=pseudotime_tf, y=simulation, color = paste0(gene_pair_2[2], " (sim)")), lwd=0.75) +
  scale_color_manual(
    name = "gene",
    values = setNames(
      c("#E64B35", "#95A3F3", "#0D1B6A"),
      c(gene_pair_2[1], gene_pair_2[2], paste0(gene_pair_2[2], " (sim)"))
    )
  ) + theme_bw() + 
  ylab('log10(counts+1)') + 
  xlab("pseudotime") + 
  guides(colour="none")



### ORIGINAL TRAJECTORIES ####
gene_class = c(rep("TF", 11), rep("target", 13))
pseudotimeDE_curves = vector(mode="list", length=length(gene_vec))
for(g in seq(gene_vec)) {
  
  gene = gene_vec[g]
  index = which(res$gene %in% gene) # find index of gene
  gene_curve = PseudotimeDE::plotCurve(gene.vec = gene, 
                                       ori.tbl = retro_tbl[[2]], 
                                       assay = 'logcounts',
                                       mat = retro_tbl[[1]], 
                                       model.fit = res$gam.fit[index])
  
  gene_header_col = ifelse(gene_class[g] == "TF", "#F89C86", "#95A3F3")
  gene_curve[["theme"]][["strip.background"]][["fill"]] = gene_header_col
  
  pseudotimeDE_curves[[g]] = gene_curve 
}
names(pseudotimeDE_curves) = gene_vec

pseudotimeDE_curves[["Fos"]] + (pseudotimeDE_curves[["Mt1"]] + labs(y=""))
pseudotimeDE_curves[["Egr1"]] + (pseudotimeDE_curves[["Cd9"]] + labs(y=""))

(pseudotimeDE_curves[["Jun"]] + (pseudotimeDE_curves[["Cxcl1"]] + labs(y=""))) /
   (pseudotimeDE_curves[["Nfia"]] + (pseudotimeDE_curves[["Rbp1"]] + labs(y="")))

(pseudotimeDE_curves[["Jund"]] + (pseudotimeDE_curves[["Mmp2"]] + labs(y=""))) / 
  (pseudotimeDE_curves[["Nfia"]] + (pseudotimeDE_curves[["Fn1"]] + labs(y="")))

(pseudotimeDE_curves[["Tgfb1"]] + (pseudotimeDE_curves[["Col1a2"]] + labs(y=""))) /
  (pseudotimeDE_curves[["Ctnnb1"]] + (pseudotimeDE_curves[["Dpep1"]] + labs(y="")))

(pseudotimeDE_curves[["Ctnnb1"]] + (pseudotimeDE_curves[["Emp1"]] + labs(y=""))) /
  (pseudotimeDE_curves[["Ctnnb1"]] + (pseudotimeDE_curves[["Wisp2"]] + labs(y="")))


(pseudotimeDE_curves[["Tcf4"]] + (pseudotimeDE_curves[["Acta2"]] + labs(y=""))) /
  (pseudotimeDE_curves[["Thra"]] + (pseudotimeDE_curves[["Col6a1"]] + labs(y="")))

(pseudotimeDE_curves[["Ar"]] + (pseudotimeDE_curves[["Lpl"]] + labs(y=""))) /
  (pseudotimeDE_curves[["Zeb1"]] + (pseudotimeDE_curves[["Tagln2"]] + labs(y="")))


###### ARC LENGTH TRAJECTORIES ######

arclength = fitting[[1]]
retro_tbl_arc = preproc_data(scd_myo5, arclength)

knots = c(70, 150, 200)
res_arc <- PseudotimeDE::runPseudotimeDE(gene.vec = pseudotimeDE_genes,
                                         ori.tbl = retro_tbl_arc[["ori"]],
                                         sub.tbl = retro_tbl_arc[["sub"]], ## To save time, use 100 subsamples
                                         mat = retro_tbl_arc[[1]], ## You can also use a matrix or SeuratObj as input
                                         mc.cores = 4,
                                         assay.use = "logcounts",
                                         model = "nb", # bc log-transformed
                                         knots=knots,
                                         k=length(knots))

pseudotimeDE_curves_arc = vector(mode="list", length=length(pseudotimeDE_genes))
for(g in seq(pseudotimeDE_genes)) {
  
  gene = pseudotimeDE_genes[g]
  index = which(res_arc$gene %in% gene) # find index of gene
  gene_curve = PseudotimeDE::plotCurve(gene.vec = gene, 
                                       ori.tbl = retro_tbl_arc[[2]], 
                                       assay = 'logcounts',
                                       mat = retro_tbl_arc[[1]], 
                                       model.fit = res_arc$gam.fit[index])
  
  gene_header_col = ifelse(gene_class[g] == "TF", "#F89C86", "#95A3F3")
  gene_curve[["theme"]][["strip.background"]][["fill"]] = gene_header_col
  gene_curve = gene_curve + xlab("") + ylab("")
  
  pseudotimeDE_curves_arc[[g]] = gene_curve
}
names(pseudotimeDE_curves_arc) = pseudotimeDE_genes


(pseudotimeDE_curves_arc[[1]] + pseudotimeDE_curves_arc[[2]]) / 
  (pseudotimeDE_curves_arc[[3]] + pseudotimeDE_curves_arc[[4]]) / 
  (pseudotimeDE_curves_arc[[5]] + pseudotimeDE_curves_arc[[6]]) 



###### ORDERING TRAJECTORIES ######
pseudotime_ordered = order(pseudotime)
retro_tbl_ord = preproc_data(scd_myo5, pseudotime_ordered)

knots  = c(1000,2000,3000)
res_ordered <- PseudotimeDE::runPseudotimeDE(gene.vec = pseudotimeDE_genes,
                                             ori.tbl = retro_tbl_ord[["ori"]],
                                             sub.tbl = retro_tbl_ord[["sub"]], ## To save time, use 100 subsamples
                                             mat = retro_tbl_ord[[1]], ## You can also use a matrix or SeuratObj as input
                                             mc.cores = 4,
                                             assay.use = "logcounts",
                                             model = "nb", # bc log-transformed
                                             knots=knots,
                                             k=length(knots))

pseudotimeDE_curves_ord = vector(mode="list", length=length(pseudotimeDE_genes))
for(g in seq(pseudotimeDE_genes)) {
  
  gene = pseudotimeDE_genes[g]
  index = which(res_ordered$gene %in% gene) # find index of gene
  gene_curve = PseudotimeDE::plotCurve(gene.vec = gene, 
                                       ori.tbl = retro_tbl_ord[[2]], 
                                       assay = 'logcounts',
                                       mat = retro_tbl_ord[[1]], 
                                       model.fit = res_ordered$gam.fit[index])
  
  gene_header_col = ifelse(gene_class[g] == "TF", "#F89C86", "#95A3F3")
  gene_curve[["theme"]][["strip.background"]][["fill"]] = gene_header_col
  gene_curve = gene_curve + xlab("") + ylab("")
  
  pseudotimeDE_curves_ord[[g]] = gene_curve
}


(pseudotimeDE_curves_ord[[1]] + pseudotimeDE_curves_ord[[2]]) / 
  (pseudotimeDE_curves_ord[[3]] + pseudotimeDE_curves_ord[[4]]) /
  (pseudotimeDE_curves_ord[[5]] + pseudotimeDE_curves_ord[[6]])


