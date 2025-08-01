## description

library(SingleCellExperiment)
library(PseudotimeDE)
library(tibble)
library(pracma) # isempty()

rm(list=ls())

#### Load dataset and functions ####
myelo_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/real/scd_myelo.rds"
myo_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/real/scd_myo.rds"
scd_myelo = readRDS(gzcon(url(myelo_url)))
scd_myo = readRDS(gzcon(url(myo_url)))
source("~/KaitlynRRStudio/RETRO-analysis/main_figure_scripts/tf_modeling_functions.R")

#### Myelopoiesis analysis ####

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
retro_pt_obj = readRDS(file=paste0(dir, "RETRO_PT_myelo.rds"))

# Turn in tibble object for PseudotimeDE analysis
preproc_data <- function(data, pseudotime, lineage=NULL) {
  
  x = data@assayData[["exprs"]]
  x_dim = data@experimentData@other[[1]]
  x_time = data@phenoData@data[["time"]]
  x_celltypes = data@phenoData@data[["cell_type"]]
  p = as.numeric(pseudotime)
  
  # Subset data according to one lineage
  if(!isempty(lineage)) {
    x <- x[,lineage]
    x_time <- x_time[lineage]
    x_dim <- x_dim[lineage,]
    x_celltypes <- x_celltypes[lineage]
    p <- as.numeric(pseudotime)[lineage] # already scaled to real time
  }
  
  # Convert eset object to SCE for PseudotimeDE
  data_sce <- SingleCellExperiment(as.matrix(x))
  data_sce@colData@listData[["cell_labels"]] <- x_celltypes
  data_sce@colData@listData[["time"]] <- x_time
  
  reducedDims(data_sce) <- list("PCA"=x_dim)
  assayNames(data_sce) <- 'logcounts'
  
  ## Number of subsamples
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


##### Arc length #####
lin_membership = retro_pt_obj@lin_membership
fitting = retro_pt_obj@arc_length
pseudotime = retro_pt_obj@pseudotime

retro_tbl_mac = preproc_data(scd_myelo, fitting[[3]], lineage=lin_membership[[3]]) # mac
retro_tbl_mon = preproc_data(scd_myelo, fitting[[4]], lineage=lin_membership[[4]]) # mon
retro_tbl_meg = preproc_data(scd_myelo, fitting[[1]], lineage=lin_membership[[1]]) # meg
retro_tbl_mast = preproc_data(scd_myelo, fitting[[2]], lineage=lin_membership[[2]]) # 

##### Ordered pseudotime #####
pseudotime_ordering = order(pseudotime)
retro_tbl_mac = preproc_data(scd_myelo, pseudotime_ordering, lineage=lin_membership[[3]]) # mac
retro_tbl_mon = preproc_data(scd_myelo, pseudotime_ordering, lineage=lin_membership[[4]]) # mon
retro_tbl_meg = preproc_data(scd_myelo, pseudotime_ordering, lineage=lin_membership[[1]]) # meg
retro_tbl_mast = preproc_data(scd_myelo, pseudotime_ordering, lineage=lin_membership[[2]]) # 

##### Redo PseudotimeDE Fitting #####
mega_gene_vec = c(TF = "GTF2I", target = "HSPA5")
mac_gene_vec = c(TF = "CEBPB", target = "ITGB2")
mon_gene_vec = c(TF = "SPI1", target="S100A9")
mast_gene_vec = c(TF = "POU5F1", target = "TDGF1")

knots = c(1000,2000,3000) # ordered pseudotime
# knots = c(5,20,35) # arc length

res_1 <- PseudotimeDE::runPseudotimeDE(gene.vec = mega_gene_vec,
                                       ori.tbl = retro_tbl_meg[["ori"]],
                                       sub.tbl = retro_tbl_meg[["sub"]], ## To save time, use 100 subsamples
                                       mat = retro_tbl_meg[[1]], ## You can also use a matrix or SeuratObj as input
                                       mc.cores = 4,
                                       assay.use = "logcounts",
                                       model = "nb", # bc log-transformed
                                       knots=knots,
                                       k=length(knots))

res_2 <- PseudotimeDE::runPseudotimeDE(gene.vec = mac_gene_vec,
                                       ori.tbl = retro_tbl_mac[["ori"]],
                                       sub.tbl = retro_tbl_mac[["sub"]], ## To save time, use 100 subsamples
                                       mat = retro_tbl_mac[[1]], ## You can also use a matrix or SeuratObj as input
                                       mc.cores = 4,
                                       assay.use = "logcounts",
                                       model = "nb", # bc log-transformed
                                       knots=knots,
                                       k=length(knots))

res_3 <- PseudotimeDE::runPseudotimeDE(gene.vec = mon_gene_vec,
                                       ori.tbl = retro_tbl_mon[["ori"]],
                                       sub.tbl = retro_tbl_mon[["sub"]], ## To save time, use 100 subsamples
                                       mat = retro_tbl_mon[[1]], ## You can also use a matrix or SeuratObj as input
                                       mc.cores = 4,
                                       assay.use = "logcounts",
                                       model = "nb", # bc log-transformed
                                       knots=knots,
                                       k=length(knots))

res_4 <- PseudotimeDE::runPseudotimeDE(gene.vec = mast_gene_vec,
                                       ori.tbl = retro_tbl_mast[["ori"]],
                                       sub.tbl = retro_tbl_mast[["sub"]], ## To save time, use 100 subsamples
                                       mat = retro_tbl_mast[[1]], ## You can also use a matrix or SeuratObj as input
                                       mc.cores = 4,
                                       assay.use = "logcounts",
                                       model = "nb", # bc log-transformed
                                       knots=knots,
                                       k=length(knots))


#### Plotting PseudotimeDE Curves ####

###### (1) Megakaryocyte Lineage ###### 
gene_class = c("TF", "target")
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
  
  curves_1[[i]] = gene_curve + labs(x="", y="")
}

curves_1
names(curves_1) = gene_class

###### (2) Macrophage Lineage ###### 
gene_class = c("TF", "target")
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
  
  curves_2[[i]] = gene_curve + labs(x="", y="")
}

curves_2
names(curves_2) = gene_class

###### (3) Monocyte Lineage ###### 
gene_class = c("TF", "target")
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
  
  curves_3[[i]] = gene_curve + labs(x="", y="")
}

curves_3
names(curves_3) = gene_class


###### (4) Mast Cell Lineage ###### 
gene_class = c("TF", "target")
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
  
  curves_4[[i]] = gene_curve + labs(x="", y="")
}

curves_4
names(curves_4) = gene_class

((curves_1[[1]] + curves_1[[2]]) |
    (curves_2[[1]] + curves_2[[2]])) /
  ((curves_3[[1]] + curves_3[[2]]) |
     (curves_4[[1]] + curves_4[[2]]))


#### Myocardial infarction analysis ####

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
retro_pt_obj = readRDS(file=paste0(dir, "RETRO_PT_myo.rds"))
fitting = unlist(retro_pt_obj@arc_length)
pseudotime = retro_pt_obj@pseudotime
pseudotimeDE_genes = c("Jun", "Cxcl1", "Nfia", "Rbp1")
gene_class = c("TF", "target", "TF", "target")

##### Arc length #####
retro_tbl_arc = preproc_data(scd_myo, fitting)

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
  (pseudotimeDE_curves_arc[[3]] + pseudotimeDE_curves_arc[[4]]) 



##### Ordering #####
pseudotime_ordering = order(pseudotime)
retro_tbl_ord = preproc_data(scd_myo, pseudotime_ordering)

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
  (pseudotimeDE_curves_ord[[3]] + pseudotimeDE_curves_ord[[4]]) 


