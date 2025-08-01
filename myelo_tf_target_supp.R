rm(list=ls())

# Load library/functions
library(SingleCellExperiment)
library(PseudotimeDE)
library(tibble)
source("~/KaitlynRRStudio/PseudotimeProject/TF_Modeling_Functions.R")

# Load PseudotimeDE object
load("~/KaitlynRRStudio/RETRO-analysis/tf_target_data/retro_tbl_mac.rda")
load("~/KaitlynRRStudio/RETRO-analysis/tf_target_data/retro_tbl_mon.rda")
load("~/KaitlynRRStudio/RETRO-analysis/tf_target_data/retro_tbl_mast.rda")
load("~/KaitlynRRStudio/RETRO-analysis/tf_target_data/retro_tbl_meg.rda")

# Load trajectories 
load("~/KaitlynRRStudio/RETRO-analysis/tf_target_data/myelo_ptde_traj_mac_lin.rda")
load("~/KaitlynRRStudio/RETRO-analysis/tf_target_data/myelo_ptde_traj_mon_lin.rda")
load("~/KaitlynRRStudio/RETRO-analysis/tf_target_data/myelo_ptde_traj_mast_lin.rda")
load("~/KaitlynRRStudio/RETRO-analysis/tf_target_data/myelo_ptde_traj_meg_lin.rda")

###### (1) Megakaryocyte Lineage ###### 
mega_gene_vec = c("FLI1", "CCND3", "PTTG1", "S100A4")
gene_class = c("TF", "target", "TF", "target")
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

###### (2) Macrophage Lineage ###### 
mac_gene_vec = c("SPI1", "TYROBP", "JUND", "CSTA")
gene_class = c("TF", "target", "TF", "target")
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

###### (3) Monocyte Lineage ###### 
mon_gene_vec = c("SPI1", "LSP1", "CEBPD", "ALOX5AP")
gene_class = c("TF", "target", "TF", "target")
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


###### (4) Mast cell Lineage ###### 
mast_gene_vec = c("GATA1", "SPI1", "ZEB2", "VIM")
gene_class = c("TF", "target", "TF", "target")
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


###### (5) Progenitor Cell State ###### 
prog_gene_vec = c("POU5F1", "TDGF1", "POU5F1", "BMP4")
gene_class = c("TF", "target", "TF", "target")
res_5 = tf_target_comp_mast[["results"]] # pseudotimeDE param fitting
index_5 = match(prog_gene_vec, res_5$gene)  # indexing

curves_5 = vector(mode="list", length=length(prog_gene_vec))
for(i in seq(index_5)) {
  
  ind = index_5[i]
  gene_curve <- PseudotimeDE::plotCurve(gene.vec = res_5$gene[ind],
                                        ori.tbl = retro_tbl_mast[[2]],
                                        assay = 'logcounts',
                                        mat = retro_tbl_mast[[1]],
                                        model.fit = res_5$gam.fit[ind])
  
  gene_header_col = ifelse(gene_class[i] == "TF", "#F89C86", "#95A3F3")
  gene_curve[["theme"]][["strip.background"]][["fill"]] = gene_header_col
  
  curves_5[[i]] = gene_curve
}

curves_5
names(curves_5) = gene_class


##### PLOTTING #####

(curves_1[[1]] + curves_1[[2]]) | (curves_1[[3]] + curves_1[[4]])  # mega
(curves_2[[1]] + curves_2[[2]]) | (curves_2[[3]] + curves_2[[4]])  # macro
(curves_3[[1]] + curves_3[[2]]) | (curves_3[[3]] + curves_3[[4]])  # monocyte
(curves_4[[1]] + curves_4[[2]]) | (curves_4[[3]] + curves_4[[4]])  # mast cell
(curves_5[[1]] + curves_5[[2]]) | (curves_5[[3]] + curves_5[[4]])  # progenitor




