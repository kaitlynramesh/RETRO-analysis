
# load datasets
cyclic_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/benchmark/RETRO_tree_v51.rda"
scd_twocycles = readRDS(gzcon(url(cyclic_url)))


# load scoring objects from cluster (all retro_obj)



retro_obj <- get_num_lineages(retro_obj, percent=0.1, cutoff=0.8, threshold=0.8) # top MST 

retro_obj <- get_bezier_curve(retro_obj, extension=2) # MST --> Bézier curve

retro_pt_obj <- get_mapped_cells(retro_obj)
retro_pt_obj <- pseudotime_fit(retro_pt_obj, retro_obj)
graph_list <- retro_pt_obj$Graph_List # e vs lambda
pseudotime <- retro_pt_obj$Pseudotime

bcurves <- sapply(retro_obj[["Curve"]], "[", 3)
# get time
# get cell type

g0 = ggplot(as.data.frame(scd_bifV2b$PCA$x), aes(x=PC1, y=PC2, colour=Time)) + 
  geom_point() + 
  guides(colour=guide_legend(ncol=1)) + 
  theme_bw()

Pseudotime = as.numeric(pseudotime)
g1 = ggplot(as.data.frame(coordinates), aes(x=PC1, y=PC2, colour=Pseudotime)) + 
  geom_point() + 
  # scale_colour_viridis_c(option="inferno") + 
  scale_colour_continuous(low="blue", high="orange") + 
  geom_point(data = as.data.frame(bcurves[[1]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(bcurves[[2]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()


dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
lin_membership = retro_pt_obj$Lin_Membership
fitting = retro_pt_obj$Fitting

save(lin_membership, file=paste0(dir, "retro_bifV2b_linmem.rda"))
save(fitting, file=paste0(dir, "retro_bifV2b_fitting.rda"))
save(pseudotime, file=paste0(dir, "retro_bifV2b_pseudotime.rda"))


