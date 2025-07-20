### SUPPLEMENTAL PLOTS ###

# Plotting all TF-target pairs from main figure
# all gene pairs
d0_d3_phase = cbind(c("Fos", "Cebpb", "Nfia", "Egr1"),
                    c("Mt1", "Crip2", "Fn1", "Col1a2"))
d5_d7_phase = cbind(c("Sox9", "Ctnnb1", "Ctnnb1", "Ctnnb1"),
                    c("Cdkn1a", "Dpep1", "Wisp2", "Emp1"))
d14_d28_phase = cbind(c("Zeb1", "Ar", "Jund", "Tcf4"),
                      c("Tagln2", "Lpl", "Timp1", "Id2"))

all_gene_pairs = c(paste(d0_d3_phase[,1], d0_d3_phase[,2], sep="_"), 
                   paste(d5_d7_phase[,1], d5_d7_phase[,2], sep="_"), 
                   paste(d14_d28_phase[,1], d14_d28_phase[,2], sep="_"))
index = which(paste(tf_summary_mat$TF, tf_summary_mat$target,sep="_") %in% all_gene_pairs)
tf_summary_mat[index,] # index in summary matrix

knots = c(5, 10, 15)
gene_vec = unique(c(d0_d3_phase, d5_d7_phase, d14_d28_phase)) # unique genes
res <- PseudotimeDE::runPseudotimeDE(gene.vec = gene_vec,
                                     ori.tbl = retro_tbl[["ori"]],
                                     sub.tbl = retro_tbl[["sub"]], ## To save time, use 100 subsamples
                                     mat = retro_tbl[[1]], ## You can also use a matrix or SeuratObj as input
                                     mc.cores = 4,
                                     assay.use = "logcounts",
                                     model = "nb", # bc log-transformed
                                     knots=knots,
                                     k=length(knots))

gene_class = c(rep("TF",4),rep("target",4), rep("TF", 2),rep("target",4),
               rep("TF",4),rep("target",4))
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
  gene_curve = gene_curve + xlab("") + ylab("log10(count+1)")
  
  pseudotimeDE_curves[[g]] = gene_curve
}
names(pseudotimeDE_curves) = gene_vec


d0_d3_phase = cbind(c("Fos", "Cebpb", "Nfia", "Egr1"),
                    c("Mt1", "Crip2", "Fn1", "Col1a2"))
d5_d7_phase = cbind(c("Sox9", "Ctnnb1", "Ctnnb1", "Ctnnb1"),
                    c("Cdkn1a", "Dpep1", "Wisp2", "Emp1"))
d14_d28_phase = cbind(c("Zeb1", "Ar", "Jund", "Tcf4"),
                      c("Tagln2", "Lpl", "Timp1", "Id2"))


(pseudotimeDE_curves[["Fos"]] + pseudotimeDE_curves[["Mt1"]]) / 
  (pseudotimeDE_curves[["Cebpb"]] + pseudotimeDE_curves[["Crip2"]]) 

  (pseudotimeDE_curves[["Nfia"]] + pseudotimeDE_curves[["Fn1"]]) / 
  (pseudotimeDE_curves[["Egr1"]] + pseudotimeDE_curves[["Col1a2"]]) 

(pseudotimeDE_curves[["Sox9"]] + pseudotimeDE_curves[["Cdkn1a"]]) / 
  (pseudotimeDE_curves[["Ctnnb1"]] + pseudotimeDE_curves[["Dpep1"]]) 

(pseudotimeDE_curves[["Ctnnb1"]] + pseudotimeDE_curves[["Wisp2"]]) /
  (pseudotimeDE_curves[["Ctnnb1"]] + pseudotimeDE_curves[["Emp1"]])

(pseudotimeDE_curves[["Zeb1"]] + pseudotimeDE_curves[["Tagln2"]]) / 
  (pseudotimeDE_curves[["Ar"]] + pseudotimeDE_curves[["Lpl"]]) 

(pseudotimeDE_curves[["Jund"]] + pseudotimeDE_curves[["Timp1"]]) / 
  (pseudotimeDE_curves[["Ar"]] + pseudotimeDE_curves[["Lpl"]]) 

(pseudotimeDE_curves[["Zeb1"]] + pseudotimeDE_curves[["Tagln2"]]) / 
  (pseudotimeDE_curves[["Tcf4"]] + pseudotimeDE_curves[["Id2"]]) 



