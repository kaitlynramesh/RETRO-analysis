# Identify TF/target gene regulatory pairs from myocardial infarction data 

library(SingleCellExperiment)
library(NetAct)
data(mDB) # from TRRUST TF/target database 

# Turn in tibble object for PseudotimeDE analysis 
preproc_data <- function(data, pseudotime) {

  # obtain data for SCE
  x = data@assayData[["exprs"]]
  x_dim = data@experimentData@other[["PCA"]][[1]]
  x_time = data@phenoData@data[["time"]]
  x_celltypes = data@phenoData@data[["cell_type"]]
  p = as.numeric(pseudotime)
  
  # Convert eset object to SCE for PseudotimeDE
  data_sce <- SingleCellExperiment(as.matrix(x))
  data_sce@colData@listData[["cell_labels"]] <- x_celltypes
  data_sce@colData@listData[["time"]] <- x_time
  
  reducedDims(data_sce) <- list(PCA=x_dim)
  assayNames(data_sce) <- 'logcounts'
  
  # Number of subsamples
  n = 20
  
  index <- mclapply(seq_len(n), function(x) {
    sample(x = c(1:dim(data_sce)[2]), size = 0.8*dim(data_sce)[2], replace = FALSE)
  })
  
  ori_tbl <- tibble(cell = colnames(data_sce), pseudotime = p)
  sub_tbl <- lapply(index, function(x) {
    sce <- data_sce[,x]
    mat <- tibble(cell = colnames(sce), pseudotime = p[x])
    return(mat)
  })
  
  pseudotime_de_table <- list("sce"=data_sce, "ori"=ori_tbl, "sub"=sub_tbl, "p"=pseudotime)
  return(pseudotime_de_table)
}


#### Load datasets and pseudotime object ####
myo_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/real/scd_myo.rds"
scd_myo = readRDS(gzcon(url(myo_url)))

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
retro_pt_obj = readRDS(file=paste0(dir, "RETRO_PT_myo.rds"))
pseudotime = retro_pt_obj@pseudotime
time = retro_pt_obj@time

#### Pre-process data for PseudotimeDE analysis ####
retro_tbl = preproc_data(scd_myo, pseudotime)

# Get tf-target pairs 
tfs = names(mDB)
targets = unique(unlist(lapply(mDB, function(x) return(x))))
sce = retro_tbl[["sce"]]
tfs = tfs[which(tfs %in% rownames(sce))]
targets = targets[which(targets %in% rownames(sce))]

x = data@assayData[["exprs"]]
x_tfs = x[which(rownames(x) %in% tfs),] # obtain matrix only w/ TFs
x_targets = x[which(rownames(x) %in% targets),] # obtain matrix only w/ targets

vf <- Seurat::FindVariableFeatures(object = x_tfs)
variable_tfs = vf[order(vf$vst.variance, decreasing=TRUE),]
top_variable_tfs = rownames(variable_tfs)[1:100] 

vf <- Seurat::FindVariableFeatures(object = x_targets)
variable_targets = vf[order(vf$vst.variance, decreasing=TRUE),]
top_variable_targets = rownames(variable_targets)[1:150] 

targets = lapply(top_variable_tfs, function(tf) { # for each tf, find corresponding targets
  x_target = which(top_variable_targets %in% mDB[[tf]])
  x_target = top_variable_targets[x_target]
  return(x_target)
})
names(targets) = top_variable_tfs

# Remove TFs that don't have variable target genes
rm_tf = which(unlist(lapply(targets, function(tf) isempty(tf))))
targets = targets[-rm_tf]
all_genes = unique(c(unlist(targets), names(targets))) # remaining TF/target pairs

#### Run PseudotimeDE to get gene expression trajectories ####
knots = c(5, 10, 15) 
results <- PseudotimeDE::runPseudotimeDE(gene.vec = all_genes,
                                         ori.tbl = retro_tbl[["ori"]],
                                         sub.tbl = retro_tbl[["sub"]], ## To save time, using 20 subsamples
                                         mat = retro_tbl[[1]], ## You can also use a matrix or SeuratObj as input
                                         mc.cores = 4,
                                         assay.use = "logcounts",
                                         model = "nb", # bc log-transformed
                                         knots=knots,
                                         k=length(knots))

curves <- PseudotimeDE::plotCurve(gene.vec = results$gene,
                                  ori.tbl = retro_tbl[[2]],
                                  assay = 'logcounts',
                                  mat = retro_tbl[[1]],
                                  model.fit = results$gam.fit,
                                  ncol=5)

tf_target_comp = list("results" = results, "curves" = curves)

# Save results
dir = "~/KaitlynRRStudio/RETRO-analysis/tf_target_data"
save(retro_tbl, file=paste0(dir, "/myo_retro_tbl.rda"))
save(tf_target_comp, file=paste0(dir, "/myo_ptde_traj.rda"))


####

load(paste0(dir, "/myo_retro_tbl.rda"))
load(paste0(dir, "/myo_ptde_traj.rda"))

#### Organize TF/target gene expression trajectories based on regulatory interaxn ####
curves = tf_target_comp$curves

# Obtain gene trajectory per gene
pt = unlist(curves[["data"]][curves[["data"]]$gene=="Fos", 2])
gene_curves <- lapply(all_genes, function(gene) {
  curve_fitting = curves[["data"]][curves$data$gene==gene, c(2,4)]
  gene_traj = as.numeric(curve_fitting$fitted)
  return(gene_traj)
})
names(gene_curves) = all_genes

# List of TF gene expression trajectories 
tf_traj <- lapply(names(targets), function(tf) {
  x = gene_curves[[tf]][order(pt)]
  return(x)
})
names(tf_traj) = names(targets)

# List of target gene expression trajectories (organized by TF)
target_traj <- lapply(targets, function(genes) {
  
  gene_traj = vector(mode="list", length=length(genes))
  for(gene in genes) {
    gene_traj[[which(genes %in% gene)]] = gene_curves[[gene]][order(pt)]
  }
  names(gene_traj) = genes
  return(gene_traj)
})



#### Calculate max correlation and lag between TF/target gene trajectories ####
num_lags = 20
max_lag = 2
max_lag_for_filtering = 3
all_target_traj_corr <- lapply(seq(targets), function(i) {
  
  tf = names(targets)[i] # TF
  target_genes = targets[[tf]] # target
  
  # initialize matrices/vectors
  p_val_mat = matrix(0, nrow=length(target_genes), ncol=num_lags)
  cor_mat = matrix(0, nrow=length(target_genes), ncol=num_lags)
  peak_vec = seq(length(target_genes))
  lag_vec = seq(length(target_genes))
  filt_vec = rep(T, length = length(target_genes))
  
  for(target in target_genes) { # each target 
    
    x = tf_traj[[tf]] # TF trajectory
    y = target_traj[[tf]][[target]] # target trajectory
    
    p_val_vec = seq(num_lags) # init p_val vector 
    cor_vec = seq(num_lags) # init correlation vector 
    lags = seq(0, max_lag, length.out=num_lags) # lag times (≤2 days)
    pt_sorted = sort(pt) # sorted pseudotime values corresponding to gene exprs.
    
    j = 1
    all_shifts = seq(length(lags)) # initialize shifts
    for(lag in lags) {
      shift = which.min(abs(pt_sorted-lag)) # how many pseudotime-points over
      all_shifts[j] = shift # store value of window shift 
      
      x_lag = x[1:(length(x) - shift)] # TF expression
      y_lag = y[(shift + 1):length(y)] # target expression
      
      spear_val = cor.test(x_lag, y_lag, method="spearman", exact=FALSE) 
      
      cor_vec[j] = as.numeric(spear_val$estimate) # correlation value
      p_val_vec[j] = as.numeric(spear_val$p.value) # p-value for correlation
      j = j+1
    }
    
    target_ind = which(target_genes %in% target) # index for gene pair (TF--target)
    
    # Check gene-gene interaxns w/ max cor that DEFAULTS to max lag
    if(which.max(abs(cor_vec)) == length(all_shifts)) { 
      max_shift = which.min(abs(pt_sorted-max_lag_for_filtering))
      x_lag = x[1:(length(x) - max_shift)] # TF expression
      y_lag = y[(max_shift + 1):length(y)] # target expression
      
      s = cor.test(x_lag, y_lag, method="spearman", exact=FALSE)
      
      s = as.numeric(s$estimate) # correlation between traj with LARGER lag
      keep_gene_pair = ifelse(max(c(s, cor_vec)) != s, T, F) # keep gene pair w/ cor that does NOT default to maximum
      filt_vec[target_ind] = keep_gene_pair
    }
    
    high_corr_shift = all_shifts[which.max(abs(cor_vec))] # which lag yields high correlation
    high_corr_lag = lags[which.max(abs(cor_vec))] # the actual delay in time
    
    x_h_lag = x[1:(length(x) - high_corr_shift)] # obtain lagged TF trajectory
    pt_lag = pt_sorted[1:(length(pt) - high_corr_shift)]  # obtained lagged pseudotime
    peak_time = pt_lag[which.max(x_h_lag)] # pseudotime at maximum of lagged TF traj
    
    p_val_mat[target_ind,] = p_val_vec # all p-values across lag
    cor_mat[target_ind,] = cor_vec # correlation across lag
    peak_vec[target_ind] = peak_time # corresponding TF peak / target
    lag_vec[target_ind] = high_corr_lag
  }
  
  rownames(cor_mat) = rownames(p_val_mat) = target_genes
  colnames(cor_mat) = colnames(p_val_mat) = lags
  
  return(list("tf" = tf, 
              "target" = as.vector(target_genes),
              "lags" = lag_vec,
              "pval" = p_val_mat,
              "cor" = cor_mat,
              "peak_vec" = peak_vec,
              "filt_vec" = filt_vec))
})
tf_names = unlist(sapply(all_target_traj_corr, "[", "tf")) 
names(all_target_traj_corr) = tf_names # label list by regulating TF


#### Formatting into matrix ####
tf_summary = lapply(all_target_traj_corr, function(x) {
  
  tf = x[["tf"]]
  target_genes = x[["target"]]
  cor_mat = as.matrix(x[["cor"]])
  p_val_mat = as.matrix(x[["pval"]])
  peak_times = x[["peak_vec"]]
  lag_times = x[["lags"]]
  
  # initialize matrix per target gene
  tf_target_mat = matrix(NA, nrow=length(target_genes), ncol=6)
  
  j = 1
  for(target in target_genes) {
    
    max_cor_ind = which.max(abs(cor_mat[j,])) # index
    
    p_val = p_val_mat[j, max_cor_ind] # TF-target p-value AT MAX CORR
    max_cor = cor_mat[j, max_cor_ind]
    peak_time = peak_times[j] # peak time of TF
    lag_time = lag_times[j] # corresponding lag between TF/target
    
    tf_target_mat[j,] = c(tf, target, p_val, max_cor, lag_time, peak_time)
    j = j+1
  }
  
  return(tf_target_mat)
})

tf_summary_mat = as.data.frame(do.call(rbind, tf_summary))
colnames(tf_summary_mat) = c("TF", "target", "pval", "cor", "lag", "peak_time")
tf_summary_mat$lag = as.numeric(tf_summary_mat$lag)
tf_summary_mat$pval = as.numeric(tf_summary_mat$pval)
tf_summary_mat$cor = as.numeric(tf_summary_mat$cor)
tf_summary_mat$peak_time = as.numeric(tf_summary_mat$peak_time)

# Remove corr pairs that move w/ lag
lag_filter = unlist(sapply(all_target_traj_corr, "[", "filt_vec"))
tf_summary_mat = tf_summary_mat[lag_filter,]


dim(tf_summary_mat)


#### Histogram to see distribution of gene pair correlation against lag ####
tf_summary_hcp_mat = tf_summary_mat[abs(tf_summary_mat$cor) > 0.8 & 
                                      tf_summary_mat$pval < 0.05,]
ggplot(data.frame(lag=tf_summary_hcp_mat$lag), aes(x=lag)) + 
  geom_histogram(alpha=1, bins=15, fill=2, col=1) + 
  labs(title = "Histogram of Lags (d0-d2)")


#### Visualize variation in correlation with changes in lag
tf_summary_hcp_mat$pair = paste(tf_summary_hcp_mat$TF,
                                tf_summary_hcp_mat$target, sep="_")
ggplot(tf_summary_hcp_mat, aes(x=lag, y=pair, colour=cor)) + 
  geom_point(size=1)

# Find correct TF-target rows... (fix this???)
all_pairs = lapply(names(all_target_traj_corr), function(x) paste0(x, "_", all_target_traj_corr[[x]][["target"]]))
all_pairs = unlist(all_pairs)
tf_target_cor_mat = do.call(rbind, lapply(all_target_traj_corr, function(x) x[["cor"]]))
hcp_index = which(all_pairs %in% tf_summary_hcp_mat$pair)

all_tf_target_hcp_cor = tf_target_cor_mat[hcp_index,]
rownames(all_tf_target_hcp_cor) = tf_summary_hcp_mat$pair
colnames(all_tf_target_hcp_cor) = round(as.numeric(colnames(all_tf_target_hcp_cor)), 3)

# Separate pairs based on sign of interaction (MEAN value)
neg_pairs = which(rowMeans(all_tf_target_hcp_cor) < -0.75)
pos_pairs = which(rowMeans(all_tf_target_hcp_cor) > 0.75)

# Color scheme
white_to_blue <- colorRampPalette(c("#0000FF", "#8080FF", "#FFFFFF"))
white_to_red <- colorRampPalette(c("#FFFFFF", "#FF8080", "#FF0000"))
neg_col <- white_to_blue(50)
pos_col = white_to_red(50)

Heatmap(all_tf_target_hcp_cor[neg_pairs,], name="cor", col=neg_col, column_title = "lag time", column_title_side = "bottom",
        cluster_columns = F, cluster_rows = T) 
Heatmap(all_tf_target_hcp_cor[pos_pairs,], name="cor", col=pos_col, column_title = "lag time", column_title_side = "bottom", 
        cluster_columns = F, cluster_rows = T, row_names_gp = grid::gpar(fontsize = 10)) 



# Determine corresponding color for lag
gene_pairs = paste(tf_summary_mat$TF, tf_summary_mat$target, sep="_")
gene_pair_list = paste(c("Fos", "Egr1", "Jun", "Cebpb", "Nfia", "Nfia", "Sox9", "Ctnnb1", "Ctnnb1", "Ctnnb1", "Tcf4", "Jund", "Zeb1", "Ar"),
                       c("Mt1", "Cyr61", "Cxcl1", "Crip2", "Fn1", "Col1a2", "Cdkn1a", "Dpep1", "Emp1", "Wisp2", "Id2", "Timp1", "Tagln2", "Lpl"), sep="_")
peak = tf_summary_mat$peak_time[match(gene_pair_list, gene_pairs)] 
lag_value = tf_summary_mat$lag[match(gene_pair_list, gene_pairs)] 
n = 1.3 - (0.3 + (1 - 0.3) * (lag_value / max(lag_value))) # rescaling to avoid black hex color (1)
lag_color = gray(n)

m = data.frame(cbind(gene_pair_list, peak, lag_value, lag_color))
m[order(m$peak),]




gene.vec = c("Sox9", "Cd9")
plot(pt, gene_curves[[gene.vec[1]]], ylim=c(0,5),
     main= tf_summary_mat$cor[which(tf_summary_mat$TF==gene.vec[1] & tf_summary_mat$target==gene.vec[2])],
     sub= tf_summary_mat$peak_time[which(tf_summary_mat$TF==gene.vec[1] & tf_summary_mat$target==gene.vec[2])])
points(pt, gene_curves[[gene.vec[2]]], col="red")
legend("topright", legend = gene.vec, fill = c("black", "red"))

tf_summary_mat[tf_summary_mat$TF=="Sox9",]

