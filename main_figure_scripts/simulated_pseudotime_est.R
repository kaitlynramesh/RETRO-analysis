# Inferring pseudotime from ensemble MST scores from all SIMULATED datasets

# load datasets
bif_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_bifurcation.rds"
tree_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_tree.rds"
cycle_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_cycle.rds"
twocycles_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_twocycles.rds"
multicycle_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_multicycles.rds"

scd_bifurcation = readRDS(gzcon(url(bif_url)))
scd_tree = readRDS(gzcon(url(tree_url)))
scd_cycle = readRDS(gzcon(url(cycle_url)))
scd_twocycles = readRDS(gzcon(url(twocycles_url)))
scd_multicycle = readRDS(gzcon(url(multicycle_url)))


#### bifurcation ####

bif_scores = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/benchmark/RETRO_scores_bifurcation.rda"
load(gzcon(url(bif_scores)))

boxplot_scoring(retro_obj)

retro_obj <- get_num_lineages(retro_obj, percent=0.05, cutoff=0.8) # top MST 
retro_obj <- get_bezier_curve(retro_obj) # spline fitting

retro_pt_obj <- get_mapped_cells(retro_obj) # update lineage information
retro_pt_obj <- pseudotime_fit(retro_pt_obj) # pseudotime estimation

bcurves <- sapply(retro_pt_obj@RETRO_Curve, "[", 2)
pca_x = retro_obj@coordinates
time = as.factor(retro_obj@time)

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=time)) + 
  geom_point() + 
  guides(colour=guide_legend(ncol=1)) + 
  theme_bw()

pseudotime = retro_pt_obj@pseudotime

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=pseudotime)) + 
  geom_point() + 
  scale_colour_continuous(low="blue", high="orange") + 
  geom_point(data = as.data.frame(bcurves[[1]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(bcurves[[2]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
saveRDS(retro_pt_obj, file=paste0(dir, "RETRO_PT_bifurcation.rds"))




#### tree ####

rm(retro_obj)
rm(retro_pt_obj)

tree_scores = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/benchmark/RETRO_scores_tree.rda"
load(gzcon(url(tree_scores)))

boxplot_scoring(retro_obj)

retro_obj <- get_num_lineages(retro_obj, percent=0.05, cutoff=0.8) # top MST 
retro_obj <- get_bezier_curve(retro_obj) # spline fitting

retro_pt_obj <- get_mapped_cells(retro_obj) # update lineage information
retro_pt_obj <- pseudotime_fit(retro_pt_obj) # pseudotime estimation

bcurves <- sapply(retro_pt_obj@RETRO_Curve, "[", 2)
pca_x = retro_obj@coordinates
time = as.factor(retro_obj@time)

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=time)) + 
  geom_point() + 
  guides(colour=guide_legend(ncol=1)) + 
  theme_bw()

pseudotime = retro_pt_obj@pseudotime

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=pseudotime)) + 
  geom_point() + 
  scale_colour_continuous(low="blue", high="orange") + 
  geom_point(data = as.data.frame(bcurves[[1]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(bcurves[[2]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(bcurves[[3]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(bcurves[[4]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
saveRDS(retro_pt_obj, file=paste0(dir, "RETRO_PT_tree.rds"))



#### cycle ####

rm(retro_obj)
rm(retro_pt_obj)

cycle_scores = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/benchmark/RETRO_scores_cycle.rda"
load(gzcon(url(cycle_scores)))

boxplot_scoring(retro_obj)

retro_obj <- get_num_lineages(retro_obj, percent=0.05, cutoff=0.8) # top MST 
retro_obj <- get_bezier_curve(retro_obj) # spline fitting

retro_pt_obj <- get_mapped_cells(retro_obj) # update lineage information
retro_pt_obj <- pseudotime_fit(retro_pt_obj) # pseudotime estimation

bcurves <- sapply(retro_pt_obj@RETRO_Curve, "[", 2)
pca_x = retro_obj@coordinates
time = as.factor(retro_obj@time)

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=time)) + 
  geom_point() + 
  guides(colour=guide_legend(ncol=1)) + 
  theme_bw()

pseudotime = retro_pt_obj@pseudotime

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=pseudotime)) + 
  geom_point() + 
  scale_colour_continuous(low="blue", high="orange") + 
  geom_point(data = as.data.frame(bcurves[[1]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
saveRDS(retro_pt_obj, file=paste0(dir, "RETRO_PT_cycle.rds"))





#### two cycles ####

rm(retro_obj)
rm(retro_pt_obj)

twocycle_scores = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/benchmark/RETRO_scores_twocycles.rda"
load(gzcon(url(twocycle_scores)))

boxplot_scoring(retro_obj)

retro_obj <- get_num_lineages(retro_obj, percent=0.05, cutoff=0.8) # top MST 
retro_obj <- get_bezier_curve(retro_obj) # spline fitting

retro_pt_obj <- get_mapped_cells(retro_obj) # update lineage information
retro_pt_obj <- pseudotime_fit(retro_pt_obj) # pseudotime estimation

bcurves <- sapply(retro_pt_obj@RETRO_Curve, "[", 2)
pca_x = retro_obj@coordinates
time = as.factor(retro_obj@time)

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=time)) + 
  geom_point() + 
  guides(colour=guide_legend(ncol=1)) + 
  theme_bw()

pseudotime = retro_pt_obj@pseudotime

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=pseudotime)) + 
  geom_point() + 
  scale_colour_continuous(low="blue", high="orange") + 
  geom_point(data = as.data.frame(bcurves[[1]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
saveRDS(retro_pt_obj, file=paste0(dir, "RETRO_PT_twocycles.rds"))





#### two cycles and bifurcation ####

rm(retro_obj)
rm(retro_pt_obj)

multicycle_scores = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/benchmark/RETRO_scores_multicycles.rda"
load(gzcon(url(multicycle_scores)))

boxplot_scoring(retro_obj)

retro_obj <- get_num_lineages(retro_obj, percent=0.05, cutoff=0.8) # top MST 
retro_obj <- get_bezier_curve(retro_obj) # spline fitting

retro_pt_obj <- get_mapped_cells(retro_obj) # update lineage information
retro_pt_obj <- pseudotime_fit(retro_pt_obj) # pseudotime estimation

bcurves <- sapply(retro_pt_obj@RETRO_Curve, "[", 2)
pca_x = retro_obj@coordinates
time = as.factor(retro_obj@time)

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=time)) + 
  geom_point() + 
  guides(colour=guide_legend(ncol=1)) + 
  theme_bw()

pseudotime = retro_pt_obj@pseudotime

ggplot(as.data.frame(pca_x), aes(x=PC1, y=PC2, colour=pseudotime)) + 
  geom_point() + 
  scale_colour_continuous(low="blue", high="orange") + 
  geom_point(data = as.data.frame(bcurves[[1]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(bcurves[[2]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
saveRDS(retro_pt_obj, file=paste0(dir, "RETRO_PT_multicycle.rds"))




