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
dir = "~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/"
load(paste0(dir, "sctda_bifV2b_pseudotime.rda")) # sctda
scTDA = sctda_list[["scTDA"]]
sctime = sctda_list[["sctime"]]
sctime = sctime[-length(sctime)] # remove NA value
scTDA = scTDA[-length(scTDA)] # remove NA value
load(paste0(dir, "slingshot_bifV2b_pseudotime.rda")) # slingshot
load(paste0(dir, "psuper_bifV2b_pseudotime.rda")) # psupertime_pt
load(paste0(dir, "retro_bifV2b_pseudotime.rda")) # retro

# Load data
load("~/KaitlynRRStudio/PseudotimeProject/scData_Scoring/scd_bifV2b.rda")
time <- scd_bifV2b$Time

# mean pseudotime over time
mean_retro <- unlist(lapply(split(pseudotime, time), function(x) mean(x)))
mean_sling <- unlist(lapply(split(slingshot_pt, time), function(x) mean(x)))
mean_psuper <- unlist(lapply(split(psupertime_pt, time), function(x) mean(x)))
mean_scTDA <- unlist(lapply(split(scTDA, sctime), function(x) mean(x)))

av_time <- c(names(mean_retro), names(mean_sling), names(mean_psuper), names(mean_scTDA))
av_time <- as.numeric(av_time)

mean_pseudotime <- data.frame(c(mean_retro, mean_sling, mean_psuper, mean_scTDA))
group <- c(rep('RETRO', length(mean_retro)),
           rep('Slingshot', length(mean_sling)), 
           rep('psupertime_pt', length(mean_psuper)), 
           rep('scTDA', length(mean_scTDA)))
mean_pseudotime <- cbind(av_time, mean_pseudotime, as.factor(group))
colnames(mean_pseudotime) <- c('Time', 'Mean_Pseudotime', 'Algorithm')

ggplot(data=mean_pseudotime, aes(x=Time, y=Mean_Pseudotime, fill=Algorithm)) +
  geom_line(aes(col=Algorithm), linewidth=1) + 
  geom_point(aes(col=Algorithm), size=2) + 
  scale_color_brewer(palette='Spectral') + 
  theme_bw() + 
  ggtitle('Average Inferred Pseudotime for Bifurcating Simulated Data')

retro_df <- split_time_data(time, pseudotime)
sling_df <- split_time_data(time, slingshot_pt)
psuper_df <- split_time_data(time, psupertime_pt)
scTDA_df <- split_time_data(sctime, scTDA)

time <- scd_bifV2b$Time
colnames(retro_df)[2] = "RETRO"
g1 <- ggplot(retro_df, aes(x=time, y=RETRO)) + 
  geom_boxplot(col='red') + geom_point(col='red') + 
  ggtitle('Real Time vs. RETRO Pseudotime') + theme_bw()

colnames(sling_df)[2] = "Slingshot"
g2 <- ggplot(sling_df, aes(x=time, y=Slingshot)) + 
  geom_boxplot(col='darkblue') + geom_point(col='darkblue') +
  ggtitle('Real Time vs. Slingshot') + theme_bw()

colnames(psuper_df)[2] = "psupertime_pt"
g3 <- ggplot(psuper_df, aes(x=time, y=psupertime_pt)) + 
  geom_boxplot(col='darkgreen') + geom_point(col='darkgreen') +
  ggtitle('Real Time vs. psupertime_pt') + theme_bw()

colnames(scTDA_df)[2] = "scTDA"
g4 <- ggplot(scTDA_df, aes(x=time, y=scTDA)) + 
  geom_boxplot(col='purple4') + geom_point(col='purple4') +
  ggtitle('Real Time vs. scTDA') + theme_bw()

ggarrange(g1,g2,g3,g4, ncol=2, nrow=2)
bif_boxplots = list(g1, g2, g3, g4)
save(bif_boxplots, file="~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/bifV2b_boxplots.rda")

# plot all boxplots
rm(list=ls())
load("~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/bifV2b_boxplots.rda")
load("~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/multicyclic_boxplots.rda")
load("~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/twocyclic_boxplots.rda")

ggarrange(bif_boxplots[[1]],  bif_boxplots[[2]], bif_boxplots[[3]], bif_boxplots[[4]], nrow=2, ncol=2)
ggarrange(twocyclic_boxplots[[1]],  twocyclic_boxplots[[2]], twocyclic_boxplots[[3]], twocyclic_boxplots[[4]], nrow=2, ncol=2)
ggarrange(multicyclic_boxplots[[1]],  multicyclic_boxplots[[2]], multicyclic_boxplots[[3]], multicyclic_boxplots[[4]], nrow=2, ncol=2)

# Stored values
load("~/KaitlynRRStudio/PseudotimeProject/scData_Scoring/scd_bifV2b.rda")
dir = "~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/"
load(paste0(dir, "retro_bifV2b_pseudotime.rda")) # retro
time = as.factor(scd_bifV2b$Time)
pca_x = scd_bifV2b$PCA$x
g_time = ggplot(as.data.frame(pca_x), aes(x=PC1,y=PC2, colour=time)) + 
  geom_point() + 
  theme_bw() +
  theme(legend.position = "right")  
# guides(colour=guide_legend(ncol=3))
g_pseudotime = ggplot(as.data.frame(pca_x), aes(x=PC1,y=PC2, colour=pseudotime)) + 
  geom_point() + 
  scale_colour_continuous(low="blue", high="orange") + 
  theme_bw() +
  theme(legend.position = "right")
ggarrange(g_time, g_pseudotime, ncol=2)

library(patchwork)

boxplots = ggarrange(twocyclic_boxplots[[1]],  twocyclic_boxplots[[2]], twocyclic_boxplots[[3]], twocyclic_boxplots[[4]], nrow=2, ncol=2)

g_time + g_pseudotime | boxplots





#### Two cycles ####

# Stored values
dir = "~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/"
load(paste0(dir, "sctda_twocyclicV3_pseudotime.rda")) # sctda
scTDA = sctda_list[["scTDA"]]
sctime = sctda_list[["sctime"]]
load(paste0(dir, "slingshot_twocyclicV3_pseudotime.rda")) # slingshot
load(paste0(dir, "psuper_twocyclicV3_pseudotime.rda")) # psupertime_pt
load(paste0(dir, "retro_twocyclicV3_pseudotime.rda")) # retro

# Load data
load("~/KaitlynRRStudio/PseudotimeProject/scData_Scoring/scd_twocyclicV3.rda")
conttime <- scd_twocyclicV3$Time_Cont
time <- scd_twocyclicV3$Time

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
corr_test <- cbind(corr_test, rep(c('RETRO', 'psupertime_pt', 'scTDA', 'Slingshot'), 2))
corr_test <- cbind(corr_test, c(rep('Kendall_Tau', 4), rep('Spearman', 4)))
colnames(corr_test) <- c('Correlation', 'Algorithm', 'Test')

ggplot(data=corr_test, aes(x=Test, y=Correlation, fill=Algorithm)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_brewer(palette='Spectral') + 
  theme_bw() + 
  ggtitle('Inferred Pseudotime Correlation with Simulated Time')


# mean pseudotime over time
mean_retro <- unlist(lapply(split(pseudotime, time), function(x) mean(x)))
mean_sling <- unlist(lapply(split(slingshot_pt, time), function(x) mean(x)))
mean_psuper <- unlist(lapply(split(psupertime_pt, time), function(x) mean(x)))
mean_scTDA <- unlist(lapply(split(scTDA, sctime), function(x) mean(x)))

av_time <- c(names(mean_retro), names(mean_sling), names(mean_psuper), names(mean_scTDA))
av_time <- as.numeric(av_time)

mean_pseudotime <- data.frame(c(mean_retro, mean_sling, mean_psuper, mean_scTDA))
group <- c(rep('RETRO', length(mean_retro)),
           rep('Slingshot', length(mean_sling)), 
           rep('psupertime_pt', length(mean_psuper)), 
           rep('scTDA', length(mean_scTDA)))
mean_pseudotime <- cbind(av_time, mean_pseudotime, as.factor(group))
colnames(mean_pseudotime) <- c('Time', 'Mean_Pseudotime', 'Algorithm')

ggplot(data=mean_pseudotime, aes(x=Time, y=Mean_Pseudotime, fill=Algorithm)) +
  geom_line(aes(col=Algorithm), linewidth=1) + 
  geom_point(aes(col=Algorithm), size=2) + 
  scale_color_brewer(palette='Spectral') + 
  theme_bw() + 
  ggtitle('Average Inferred Pseudotime for Cyclic Simulated Data')



retro_df <- split_time_data(time, pseudotime)
sling_df <- split_time_data(time, slingshot_pt)
psuper_df <- split_time_data(time, psupertime_pt)
scTDA_df <- split_time_data(sctime, scTDA)

time <- scd_twocyclicV3$Time
colnames(retro_df)[2] = "RETRO"
g1 <- ggplot(retro_df, aes(x=time, y=RETRO)) + 
  geom_boxplot(col='red') + geom_point(col='red') + 
  ggtitle('Real Time vs. RETRO Pseudotime') + theme_bw()

colnames(sling_df)[2] = "Slingshot"
g2 <- ggplot(sling_df, aes(x=time, y=Slingshot)) + 
  geom_boxplot(col='darkblue') + geom_point(col='darkblue') +
  ggtitle('Real Time vs. Slingshot') + theme_bw()

colnames(psuper_df)[2] = "psupertime_pt"
g3 <- ggplot(psuper_df, aes(x=time, y=psupertime_pt)) + 
  geom_boxplot(col='darkgreen') + geom_point(col='darkgreen') +
  ggtitle('Real Time vs. psupertime_pt') + theme_bw()

colnames(scTDA_df)[2] = "scTDA"
g4 <- ggplot(scTDA_df, aes(x=time, y=scTDA)) + 
  geom_boxplot(col='purple4') + geom_point(col='purple4') +
  ggtitle('Real Time vs. scTDA') + theme_bw()

ggarrange(g1,g2,g3,g4, ncol=2, nrow=2)
twocyclic_boxplots = list(g1, g2, g3, g4)
save(twocyclic_boxplots, file="~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/twocyclic_boxplots.rda")





#### Two cycles and bifurcation #### 

# Stored values
dir = "~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/"
load(paste0(dir, "sctda_multicyclic2b_pseudotime.rda")) # sctda
scTDA = sctda_list[["scTDA"]]
sctime = sctda_list[["sctime"]]
sctime = sctime[-length(sctime)] # remove final value (NA)
scTDA = scTDA[-length(scTDA)]
load(paste0(dir, "slingshot_multicyclic2b_pseudotime.rda")) # slingshot
load(paste0(dir, "psuper_multicyclic2b_pseudotime.rda")) # psupertime_pt
load(paste0(dir, "retro_multicyclic2b_pseudotime.rda")) # retro

# Load data
load("~/KaitlynRRStudio/PseudotimeProject/scData_Scoring/scd_multicyclic2b.rda")
conttime <- scd_multicyclic2b$Time_Cont
time <- scd_multicyclic2b$Time

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


# mean pseudotime over time
mean_retro <- unlist(lapply(split(pseudotime, time), function(x) mean(x)))
mean_sling <- unlist(lapply(split(slingshot_pt, time), function(x) mean(x)))
mean_psuper <- unlist(lapply(split(psupertime_pt, time), function(x) mean(x)))
mean_scTDA <- unlist(lapply(split(scTDA, sctime), function(x) mean(x)))

av_time <- c(names(mean_retro), names(mean_sling), names(mean_psuper), names(mean_scTDA))
av_time <- as.numeric(av_time)

mean_pseudotime <- data.frame(c(mean_retro, mean_sling, mean_psuper, mean_scTDA))
group <- c(rep('RETRO', length(mean_retro)),
           rep('Slingshot', length(mean_sling)), 
           rep('psupertime_pt', length(mean_psuper)), 
           rep('scTDA', length(mean_scTDA)))
mean_pseudotime <- cbind(av_time, mean_pseudotime, as.factor(group))
colnames(mean_pseudotime) <- c('Time', 'Mean_Pseudotime', 'Algorithm')

ggplot(data=mean_pseudotime, aes(x=Time, y=Mean_Pseudotime, fill=Algorithm)) +
  geom_line(aes(col=Algorithm), linewidth=1) + 
  geom_point(aes(col=Algorithm), size=2) + 
  scale_color_brewer(palette='Spectral') + 
  theme_bw() + 
  ggtitle('Average Inferred Pseudotime for Cyclic Simulated Data')



retro_df <- split_time_data(time, pseudotime)
sling_df <- split_time_data(time, slingshot_pt)
psuper_df <- split_time_data(time, psupertime_pt)
scTDA_df <- split_time_data(sctime, scTDA)

time <- scd_multicyclic2b$Time
colnames(retro_df)[2] = "RETRO"
g1b <- ggplot(retro_df, aes(x=time, y=RETRO)) + 
  geom_boxplot(col='red') + geom_point(col='red') + 
  ggtitle('Real Time vs. RETRO Pseudotime') + theme_bw()

colnames(sling_df)[2] = "Slingshot"
g2b <- ggplot(sling_df, aes(x=time, y=Slingshot)) + 
  geom_boxplot(col='darkblue') + geom_point(col='darkblue') +
  ggtitle('Real Time vs. Slingshot') + theme_bw()

colnames(psuper_df)[2] = "psupertime_pt"
g3b <- ggplot(psuper_df, aes(x=time, y=psupertime_pt)) + 
  geom_boxplot(col='darkgreen') + geom_point(col='darkgreen') +
  ggtitle('Real Time vs. psupertime_pt') + theme_bw()

colnames(scTDA_df)[2] = "scTDA"
g4b <- ggplot(scTDA_df, aes(x=time, y=scTDA)) + 
  geom_boxplot(col='purple4') + geom_point(col='purple4') +
  ggtitle('Real Time vs. scTDA') + theme_bw()

ggarrange(g1b,g2b,g3b,g4b,ncol=2, nrow=2)
multicyclic_boxplots = list(g1b, g2b, g3b, g4b)
save(multicyclic_boxplots, file="~/KaitlynRRStudio/PseudotimeProject/Pseudotime_Benchmark/multicyclic_boxplots.rda")



