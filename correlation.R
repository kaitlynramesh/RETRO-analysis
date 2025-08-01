# Correlation table

#### RETRO analysis ####
scd_bifurcation = readRDS("~/KaitlynRRStudio/RETRO-analysis/synthetic/scd_bifurcation.rds")
scd_tree = readRDS("~/KaitlynRRStudio/RETRO-analysis/synthetic/scd_tree.rds")
scd_cycle = readRDS("~/KaitlynRRStudio/RETRO-analysis/synthetic/scd_cycle.rds")

RETRO_PT_bifurcation <- readRDS("~/KaitlynRRStudio/RETRO-analysis/benchmark/RETRO_PT_bifurcation.rds")
RETRO_PT_cycle <- readRDS("~/KaitlynRRStudio/RETRO-analysis/benchmark/RETRO_PT_cycle.rds")
RETRO_PT_tree <- readRDS("~/KaitlynRRStudio/RETRO-analysis/benchmark/RETRO_PT_tree.rds")
RETRO_PT_twocycles <- readRDS("~/KaitlynRRStudio/RETRO-analysis/benchmark/RETRO_PT_twocycles.rds")
RETRO_PT_multicycle <- readRDS("~/KaitlynRRStudio/RETRO-analysis/benchmark/RETRO_PT_multicycle.rds")

pseudotime = RETRO_PT_bifurcation@pseudotime
cont_time = scd_bifurcation@phenoData@data[["cont_time"]]

cor.test(pseudotime, cont_time, method="kendall")
cor.test(pseudotime, cont_time, method="spearman")


pseudotime = RETRO_PT_cycle@pseudotime
cont_time = scd_cycle@phenoData@data[["cont_time"]]

cor.test(pseudotime, cont_time, method="kendall")
cor.test(pseudotime, cont_time, method="spearman")


pseudotime = RETRO_PT_tree@pseudotime
cont_time = scd_tree@phenoData@data[["cont_time"]]

cor.test(pseudotime, cont_time, method="kendall")
cor.test(pseudotime, cont_time, method="spearman")



pseudotime = RETRO_PT_twocycles@pseudotime
time = RETRO_PT_twocycles@time

cor.test(pseudotime, time, method="kendall")
cor.test(pseudotime, time, method="spearman")


pseudotime = RETRO_PT_multicycle@pseudotime
time = RETRO_PT_multicycle@time

cor.test(pseudotime, time, method="kendall")
cor.test(pseudotime, time, method="spearman")


pseudotime = RETRO_PT_myelo@pseudotime
time = RETRO_PT_myelo@time

cor.test(pseudotime, time, method="kendall")
cor.test(pseudotime, time, method="spearman")


pseudotime = RETRO_PT_myo@pseudotime
time = RETRO_PT_myo@time

cor.test(pseudotime, time, method="kendall")
cor.test(pseudotime, time, method="spearman")


#### Slingshot analysis #### 







