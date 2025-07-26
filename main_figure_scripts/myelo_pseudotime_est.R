# Inferring pseudotime from ensemble MST scores from myelopoiesis dataset
# FIGURES (??)

library(RETRO)
library(ggplot2)

# load datasets
myelo_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/real/scd_myelo.rds"
scd_myelo = readRDS(gzcon(url(myelo_url)))

# load scoring information
myelo_scores = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/benchmark/RETRO_scores_myelo.rda"
load(gzcon(url(myelo_scores)))

retro_obj <- get_num_lineages(retro_obj, percent=0.05, cutoff=0.5) # top MST 
retro_obj <- get_bezier_curve(retro_obj, extension=0) # spline fitting

retro_pt_obj <- get_mapped_cells(retro_obj) # update lineage information
retro_pt_obj <- pseudotime_fit(retro_pt_obj) # pseudotime estimation

bcurves <- sapply(retro_pt_obj@RETRO_Curve, "[", 2)
umap = retro_obj@coordinates
time = retro_obj@time

ggplot(as.data.frame(umap), aes(x=UMAP2L_1, y=UMAP2L_2, colour=as.factor(time))) + 
  geom_point() + 
  guides(colour=guide_legend(ncol=1)) + 
  theme_bw() + labs(colour="time")

pseudotime = retro_pt_obj@pseudotime

ggplot(as.data.frame(umap), aes(x=UMAP2L_1, y=UMAP2L_2, colour=pseudotime)) + 
  geom_point() + 
  scale_colour_continuous(low="blue", high="orange") + 
  geom_point(data = as.data.frame(bcurves[[1]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(bcurves[[2]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(bcurves[[3]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  geom_point(data = as.data.frame(bcurves[[4]][,1:2]), aes(x = V1, y = V2), size=0.5, colour='black') +
  theme_bw()

pseudotime_density(time, pseudotime)

dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
saveRDS(retro_pt_obj, file=paste0(dir, "RETRO_PT_myelo.rds"))


