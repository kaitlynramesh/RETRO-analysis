# Boxplots and correlation
# Figure 2

library(ggplot2)
library(ggpubr)
library(viridis)
library(RColorBrewer)

split_time_data <- function(time, pseudotime) {
  v <- as.data.frame(cbind(time, pseudotime))
  v$time <- as.factor(v$time)
  return(v)
}


#### Bifurcation ####

# Stored values
dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark"
retro_pt_obj <- readRDS(paste0(dir, "/RETRO_PT_bifurcation.rds"))
time = retro_pt_obj@time
pseudotime = retro_pt_obj@pseudotime
load(paste0(dir, "/existing_methods/psuper_bifurcation_pseudotime.rda"))
load(paste0(dir, "/existing_methods/slingshot_bifurcation_pseudotime.rda"))
load(paste0(dir, "/existing_methods/scTDA_bifurcation_pseudotime.rda"))
scTDA = sctda_list[["scTDA"]]
sctime = sctda_list[["sctime"]]
sctime = sctime[-length(sctime)] # remove NA value
scTDA = scTDA[-length(scTDA)] # remove NA value

retro_df <- split_time_data(time, pseudotime)
sling_df <- split_time_data(time, slingshot_pt)
psuper_df <- split_time_data(time, psupertime_pt)
scTDA_df <- split_time_data(sctime, scTDA)

colnames(retro_df)[2] = "RETRO"
g1 <- ggplot(retro_df, aes(x=time, y=RETRO)) + 
  geom_boxplot(col='red') + geom_point(col='red') + 
  ggtitle('Real Time vs. RETRO') + theme_bw()

colnames(sling_df)[2] = "Slingshot"
g2 <- ggplot(sling_df, aes(x=time, y=Slingshot)) + 
  geom_boxplot(col='darkblue') + geom_point(col='darkblue') +
  ggtitle('Real Time vs. Slingshot') + theme_bw()

colnames(psuper_df)[2] = "Psupertime"
g3 <- ggplot(psuper_df, aes(x=time, y=Psupertime)) + 
  geom_boxplot(col='darkgreen') + geom_point(col='darkgreen') +
  ggtitle('Real Time vs. Psupertime') + theme_bw()

colnames(scTDA_df)[2] = "scTDA"
g4 <- ggplot(scTDA_df, aes(x=time, y=scTDA)) + 
  geom_boxplot(col='purple4') + geom_point(col='purple4') +
  ggtitle('Real Time vs. scTDA') + theme_bw()

ggarrange(g1,g2,g3,g4, ncol=2, nrow=2)
bif_boxplots = list(g1, g2, g3, g4)




#### Two cycles ####

# Stored values
dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark"
retro_pt_obj <- readRDS(paste0(dir, "/RETRO_PT_twocycles.rds"))
time = retro_pt_obj@time
pseudotime = retro_pt_obj@pseudotime
load(paste0(dir, "/existing_methods/psuper_twocycles_pseudotime.rda"))
load(paste0(dir, "/existing_methods/slingshot_twocycles_pseudotime.rda"))
load(paste0(dir, "/existing_methods/scTDA_twocycles_pseudotime.rda"))
scTDA = sctda_list[["scTDA"]]
sctime = sctda_list[["sctime"]]
sctime = sctime[-length(sctime)] # remove NA value
scTDA = scTDA[-length(scTDA)] # remove NA value

# correlation plot
k1 <- cor.test(time, pseudotime, method='kendall')$estimate
k2 <- cor.test(time, psupertime_pt, method='kendall')$estimate
k3 <- cor.test(sctime, scTDA, method='kendall')$estimate
k4 <- cor.test(time, slingshot_pt, method='kendall')$estimate

s1 <- cor.test(conttime, pseudotime, method='spearman', exact=FALSE)$estimate
s2 <- cor.test(conttime, psupertime_pt, method='spearman', exact=FALSE)$estimate
s3 <- cor.test(sctime, scTDA, method='spearman', exact=FALSE)$estimate
s4 <- cor.test(conttime, slingshot_pt, method='spearman', exact=FALSE)$estimate

corr_test <- as.data.frame((c(k1,k2,k3,k4, s1,s2,s3,s4)))
corr_test <- cbind(corr_test, rep(c('RETRO', 'Psupertime', 'scTDA', 'Slingshot'), 2))
corr_test <- cbind(corr_test, c(rep('Kendall_Tau', 4), rep('Spearman', 4)))
colnames(corr_test) <- c('Correlation', 'Algorithm', 'Test')

ggplot(data=corr_test, aes(x=Test, y=Correlation, fill=Algorithm)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette='Spectral') + 
  theme_bw() + 
  ggtitle('Inferred Pseudotime Correlation with Simulated Time')


retro_df <- split_time_data(time, pseudotime)
sling_df <- split_time_data(time, slingshot_pt)
psuper_df <- split_time_data(time, psupertime_pt)
scTDA_df <- split_time_data(sctime, scTDA)

colnames(retro_df)[2] = "RETRO"
g1 <- ggplot(retro_df, aes(x=time, y=RETRO)) + 
  geom_boxplot(col='red') + geom_point(col='red') + 
  ggtitle('Real Time vs. RETRO') + theme_bw()

colnames(sling_df)[2] = "Slingshot"
g2 <- ggplot(sling_df, aes(x=time, y=Slingshot)) + 
  geom_boxplot(col='darkblue') + geom_point(col='darkblue') +
  ggtitle('Real Time vs. Slingshot') + theme_bw()

colnames(psuper_df)[2] = "Psupertime"
g3 <- ggplot(psuper_df, aes(x=time, y=Psupertime)) + 
  geom_boxplot(col='darkgreen') + geom_point(col='darkgreen') +
  ggtitle('Real Time vs. Psupertime') + theme_bw()

colnames(scTDA_df)[2] = "scTDA"
g4 <- ggplot(scTDA_df, aes(x=time, y=scTDA)) + 
  geom_boxplot(col='purple4') + geom_point(col='purple4') +
  ggtitle('Real Time vs. scTDA') + theme_bw()

ggarrange(g1,g2,g3,g4, ncol=2, nrow=2)
twocyclic_boxplots = list(g1, g2, g3, g4)
save(twocyclic_boxplots, file="~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/twocyclic_boxplots.rda")





#### Two cycles and bifurcation #### 

# Stored values
dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark"
retro_pt_obj <- readRDS(paste0(dir, "/RETRO_PT_multicycle.rds"))
time = retro_pt_obj@time
pseudotime = retro_pt_obj@pseudotime
load(paste0(dir, "/existing_methods/psuper_multicycles_pseudotime.rda"))
load(paste0(dir, "/existing_methods/slingshot_multicycles_pseudotime.rda"))
load(paste0(dir, "/existing_methods/scTDA_multicycles_pseudotime.rda"))
scTDA = sctda_list[["scTDA"]]
sctime = sctda_list[["sctime"]]
sctime = sctime[-length(sctime)] # remove NA value
scTDA = scTDA[-length(scTDA)] # remove NA value

# correlation plot
k1 <- cor.test(time, pseudotime, method='kendall')$estimate
k2 <- cor.test(time, psupertime_pt, method='kendall')$estimate
k3 <- cor.test(sctime, scTDA, method='kendall')$estimate
k4 <- cor.test(time, slingshot_pt, method='kendall')$estimate

s1 <- cor.test(time, pseudotime, method='spearman', exact=FALSE)$estimate
s2 <- cor.test(time, psupertime_pt, method='spearman', exact=FALSE)$estimate
s3 <- cor.test(sctime, scTDA, method='spearman', exact=FALSE)$estimate
s4 <- cor.test(time, slingshot_pt, method='spearman', exact=FALSE)$estimate

corr_test <- as.data.frame((c(k1,k2,k3,k4, s1, s2, s3, s4)))
corr_test <- cbind(corr_test, rep(c('RETRO', 'psupertime_pt', 'scTDA', 'Slingshot'), 2))
corr_test <- cbind(corr_test, c(rep('Kendall_Tau', 4), rep('Spearman', 4)))
colnames(corr_test) <- c('Correlation', 'Algorithm', 'Test')

ggplot(data=corr_test, aes(x=Test, y=Correlation, fill=Algorithm)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette='Spectral') + 
  theme_bw() + 
  ggtitle('Inferred Pseudotime Correlation with Simulated Time')


retro_df <- split_time_data(time, pseudotime)
sling_df <- split_time_data(time, slingshot_pt)
psuper_df <- split_time_data(time, psupertime_pt)
scTDA_df <- split_time_data(sctime, scTDA)

colnames(retro_df)[2] = "RETRO"
g1b <- ggplot(retro_df, aes(x=time, y=RETRO)) + 
  geom_boxplot(col='red') + geom_point(col='red') + 
  ggtitle('Real Time vs. RETRO') + theme_bw()

colnames(sling_df)[2] = "Slingshot"
g2b <- ggplot(sling_df, aes(x=time, y=Slingshot)) + 
  geom_boxplot(col='darkblue') + geom_point(col='darkblue') +
  ggtitle('Real Time vs. Slingshot') + theme_bw()

colnames(psuper_df)[2] = "Psupertime"
g3b <- ggplot(psuper_df, aes(x=time, y=Psupertime)) + 
  geom_boxplot(col='darkgreen') + geom_point(col='darkgreen') +
  ggtitle('Real Time vs. Psupertime') + theme_bw()

colnames(scTDA_df)[2] = "scTDA"
g4b <- ggplot(scTDA_df, aes(x=time, y=scTDA)) + 
  geom_boxplot(col='purple4') + geom_point(col='purple4') +
  ggtitle('Real Time vs. scTDA') + theme_bw()

ggarrange(g1b,g2b,g3b,g4b,ncol=2, nrow=2)




# plot all boxplots
rm(list=ls())
load("~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/bifV2b_boxplots.rda")
load("~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/multicyclic_boxplots.rda")
load("~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/twocyclic_boxplots.rda")

ggarrange(bif_boxplots[[1]],  bif_boxplots[[2]], bif_boxplots[[3]], bif_boxplots[[4]], nrow=2, ncol=2)
ggarrange(twocyclic_boxplots[[1]],  twocyclic_boxplots[[2]], twocyclic_boxplots[[3]], twocyclic_boxplots[[4]], nrow=2, ncol=2)
ggarrange(multicyclic_boxplots[[1]],  multicyclic_boxplots[[2]], multicyclic_boxplots[[3]], multicyclic_boxplots[[4]], nrow=2, ncol=2)

