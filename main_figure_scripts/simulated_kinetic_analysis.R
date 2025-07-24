
library(viridis)
library(ggplot2)
library(pracma)

# Function to infer alpha (kinetic value) 
calculate_alpha <- function(coord, pseudotime, time, radius, cyclic=FALSE) {
  
  m <- as.matrix(dist(coord)) # gene exprs. dist matrix
  
  t = unique(time)[-c(1,length(unique(time)))] 
  breaks = sort(unique(c(min(pseudotime), t, max(pseudotime))))
  
  bins = cut(pseudotime, breaks=breaks, labels=seq(length(breaks)-1))
  bins = as.numeric(bins)
  
  coord_neighbors = vector(mode="list", length=nrow(coord))
  
  model <- lapply(1:nrow(coord), function(i) { 
    # Identifying neighbors for cell-i
    index <- which(m[,i] < radius, arr.ind=T) 
    
    # cyclic check!!
    if(cyclic==TRUE) {
      r = which(abs(bins[index] - bins[i]) > 2)
      if(!isempty(r)) {
        index = index[-r]
      }
    }
    neighbors <- coord[index,]
    
    X <- coord[i,] - neighbors # gene exprs differences
    Xt <- t(X) # transposed
    xtx <- Xt %*% X # product of transposed and neighbors
    Xinv <- tryCatch({
      solve(xtx)
    }, error = function(e) {
      print('non-invertible')
    })
    
    Xprod <- Xinv %*% Xt
    alpha <- Xprod %*% (pseudotime[i]-pseudotime[index]) # time differences
    alpha <- 1/norm(alpha, type='2')
    return(list("alpha" = alpha, "neighbors" = index))
  })
  
  alpha = sapply(model, "[[", 1)
  coord_neighbors = lapply(model, "[[", 2)
  
  # Determined smoothed alpha values for each neighborhood
  alpha_s <- lapply(1:length(alpha), function(i) { 
    n = unlist(coord_neighbors[[i]])
    mean_alpha <- mean(alpha[n])
    return(mean_alpha)
  })
  
  alpha_s <- log(unlist(alpha_s))
  return(alpha_s)
}

scale_time <- function(time, pseudotime) {
  or <- as.vector(range(pseudotime)) # original range
  ar <- as.vector(range(time)) # actual range
  scaled_pt <- lapply(pseudotime, function(p) ((((p-or[1])*(ar[2]-ar[1])) / (or[2]-or[1])) + ar[1]) )
  scaled_pt <- unlist(scaled_pt)
  return(scaled_pt)
}



##### Load data #####

bif_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_bifurcation.rds"
twocycles_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_cycle.rds"
multicycle_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_multicycles.rds"

scd_bifurcation = readRDS(gzcon(url(bif_url)))
scd_twocycles = readRDS(gzcon(url(twocycles_url)))
scd_multicycle = readRDS(gzcon(url(multicycle_url)))



##### Bifurcation example #####

# load RETRO pseudotime object
retro_pt_obj = readRDS("~/KaitlynRRStudio/RETRO-analysis/benchmark/RETRO_PT_bifurcation.rds")
coord = scd_bifurcation@experimentData@other[["PCA"]][["x"]][,1:2] # PCA values
time = retro_pt_obj@time # sampling time
pseudotime = retro_pt_obj@pseudotime # inferred pseudotime

# load benchmark pseudotime
dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/existing_methods/"
load(paste0(dir, "slingshot_bifurcation_pseudotime.rda"))
load(paste0(dir, "psuper_bifurcation_pseudotime.rda"))
load(paste0(dir, "scTDA_bifurcation_pseudotime.rda"))

nodeCoord = sctda_list$nodeCoord # 
nodeCoord = nodeCoord[-nrow(nodeCoord),] # remove NA coordinate
scTDA_pt = sctda_list$scTDA[-length(sctda_list$scTDA)]
sctime = sctda_list$sctime # sampling time at each cell 

r = 10

alpha_r <- calculate_alpha(coord, pseudotime, time, radius=r)
alpha_s <- calculate_alpha(coord, slingshot_pt, time, radius=r)
alpha_p <- calculate_alpha(coord, psupertime_pt, time, radius=r)
alpha_sc <- calculate_alpha(nodeCoord, scTDA_pt, sctime, radius=r)

a <- min(alpha_r, alpha_p, alpha_s, alpha_sc)
b <- max(alpha_r, alpha_p, alpha_s, alpha_sc)
y_lim = c(a, b)

# plotting alpha
m1 <- ggplot(as.data.frame(coord), aes(x = PC1, y = PC2, color = alpha_r)) + 
  geom_point() +
  scale_colour_viridis(option='H') +
  theme_bw() 

m2 <- ggplot(as.data.frame(coord), aes(x = PC1, y = PC2, color = alpha_s)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=c(a,b)) +
  theme_bw() 

m3 <- ggplot(as.data.frame(coord), aes(x = PC1, y = PC2, color = alpha_p)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=c(a,b)) +
  theme_bw() 

m4 <- ggplot(as.data.frame(nodeCoord), aes(x = PC1, y = PC2, color = alpha_sc)) +
  geom_point() +
  scale_colour_viridis(option='H', limits=c(a,b)) +
  theme_bw()

ggarrange(m1,m2,m3,m4, ncol=4, common.legend = T, legend="right")


## arc-length/time fitting plot
load(paste0(dir, "retro_bifV2b_linmem.rda"))
load(paste0(dir, "retro_bifV2b_fitting.rda"))
f = rowMeans(do.call(cbind, fitting), na.rm=TRUE) # arc length values

shared = intersect(lin_membership[[1]], lin_membership[[2]])
arclen_shared = f[shared] 
arclen_1 = f[lin_membership[[1]]] 
arclen_2 = f[lin_membership[[2]]]

fitting_shared_df = cbind(pseudotime[shared], arclen_shared)
fitting_1_df = cbind(pseudotime[lin_membership[[1]]], arclen_1)
fitting_2_df = cbind(pseudotime[lin_membership[[2]]], arclen_2)
fitting_shared_df = as.data.frame(fitting_shared_df)
fitting_1_df = as.data.frame(fitting_1_df)
fitting_2_df = as.data.frame(fitting_2_df)
colnames(fitting_1_df) = colnames(fitting_2_df) = colnames(fitting_shared_df) = 
  c("pseudotime", "arc_length")

ggplot() + 
  geom_point(data = fitting_1_df, aes(y=pseudotime, x=arc_length, colour="lin-1")) + 
  geom_point(data = fitting_2_df, aes(y=pseudotime, x=arc_length, colour="lin-2")) + 
  geom_point(data = fitting_shared_df, aes(y=pseudotime, x=arc_length, colour="shared")) +
  theme_bw() + 
  labs(colour="membership") +
  xlab("arc length") + 
  theme(legend.position = "bottom")



##### Two cycles example #####

rm(retro_pt_obj)
rm(coord)
rm(time)

# load RETRO pseudotime object
retro_pt_obj = readRDS("~/KaitlynRRStudio/RETRO-analysis/benchmark/RETRO_PT_twocycles.rds")
coord = scd_twocycles@experimentData@other[["PCA"]][["x"]][,1:2] # PCA values
time = retro_pt_obj@time # sampling time
pseudotime = retro_pt_obj@pseudotime # inferred pseudotime

# load benchmark pseudotime
dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/existing_methods/"
load(paste0(dir, "slingshot_twocycles_pseudotime.rda"))
load(paste0(dir, "psuper_twocycles_pseudotime.rda"))
load(paste0(dir, "scTDA_twocycles_pseudotime.rda"))

nodeCoord = sctda_list$nodeCoord 
nodeCoord = nodeCoord[-nrow(nodeCoord),] # remove NA coordinate
scTDA_pt = sctda_list$scTDA[-length(sctda_list$scTDA)]
sctime = sctda_list$sctime # sampling time at each cell 

r = 20

alpha_r <- calculate_alpha(coord, pseudotime, time, radius=r)
alpha_s <- calculate_alpha(coord, slingshot_pt, time, radius=r)
alpha_p <- calculate_alpha(coord, psupertime_pt, time, radius=r)
alpha_sc <- calculate_alpha(nodeCoord, scTDA_pt, sctime, radius=r)

a <- min(alpha_r, alpha_p, alpha_s, alpha_sc)
b <- max(alpha_r, alpha_p, alpha_s, alpha_sc)
y_lim = c(a, b)

# plotting alpha
m1 <- ggplot(as.data.frame(coord), aes(x = PC1, y = PC2, color = alpha_r)) + 
  geom_point() +
  scale_colour_viridis(option='H') +
  theme_bw() 

m2 <- ggplot(as.data.frame(coord), aes(x = PC1, y = PC2, color = alpha_s)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=c(a,b)) +
  theme_bw() 

m3 <- ggplot(as.data.frame(coord), aes(x = PC1, y = PC2, color = alpha_p)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=c(a,b)) +
  theme_bw() 

m4 <- ggplot(as.data.frame(nodeCoord), aes(x = PC1, y = PC2, color = alpha_sc)) +
  geom_point() +
  scale_colour_viridis(option='H', limits=c(a,b)) +
  theme_bw()

ggarrange(m1,m2,m3,m4, ncol=4, common.legend = T, legend="right")


## arc-length/time fitting plot
load(paste0(dir, "retro_twocyclicV3_fitting.rda"))
f = unlist(fitting)

pseudotime_cycle_1 = pseudotime[time <= 9]
pseudotime_cycle_2 = pseudotime[time > 9]

arclen_1 = f[time <= 9] / max(f[time <= 9]) 
arclen_2 = (f[time > 9]-min(f[time>9])) / max(f[time > 9]-min(f[time>9]))

fitting_cycle_1_df = as.data.frame(cbind(pseudotime_cycle_1, arclen_1))
fitting_cycle_2_df = as.data.frame(cbind(pseudotime_cycle_2, arclen_2))
colnames(fitting_cycle_1_df) = colnames(fitting_cycle_2_df) = c("pseudotime", "arc_length")

ggplot() + 
  geom_point(data = fitting_cycle_1_df, aes(y=pseudotime, x=arc_length, colour="cycle-1")) + 
  geom_point(data = fitting_cycle_2_df, aes(y=pseudotime, x=arc_length, colour="cycle-2")) + 
  theme_bw() + 
  labs(colour="membership") +
  xlab("normalized arc length") + 
  theme(legend.position = "bottom")





##### Multicyclic2 example #####

rm(retro_pt_obj)
rm(coord)
rm(time)

# load RETRO pseudotime object
retro_pt_obj = readRDS("~/KaitlynRRStudio/RETRO-analysis/benchmark/RETRO_PT_multicycle.rds")
coord = scd_multicycle@experimentData@other[["PCA"]][["x"]][,1:2] # PCA values
time = retro_pt_obj@time # sampling time
pseudotime = retro_pt_obj@pseudotime # inferred pseudotime

# load benchmark pseudotime
dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/existing_methods/"
load(paste0(dir, "slingshot_multicycles_pseudotime.rda"))
load(paste0(dir, "psuper_multicycles_pseudotime.rda"))
load(paste0(dir, "scTDA_multicycles_pseudotime.rda"))

nodeCoord = sctda_list$nodeCoord # 
nodeCoord = nodeCoord[-nrow(nodeCoord),] # remove NA coordinate
scTDA_pt = sctda_list$scTDA[-length(sctda_list$scTDA)]
sctime = sctda_list$sctime # sampling time at each cell 

r = 20

alpha_r <- calculate_alpha(coord, pseudotime, time, radius=r)
alpha_s <- calculate_alpha(coord, slingshot_pt, time, radius=r)
alpha_p <- calculate_alpha(coord, psupertime_pt, time, radius=r)
alpha_sc <- calculate_alpha(nodeCoord, scTDA_pt, sctime, radius=r)

a <- min(alpha_r, alpha_p, alpha_s, alpha_sc)
b <- max(alpha_r, alpha_p, alpha_s, alpha_sc)
y_lim = c(a, b)

# plotting alpha
m1 <- ggplot(as.data.frame(coord), aes(x = PC1, y = PC2, color = alpha_r)) + 
  geom_point() +
  scale_colour_viridis(option='H') +
  theme_bw() 

m2 <- ggplot(as.data.frame(coord), aes(x = PC1, y = PC2, color = alpha_s)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=c(a,b)) +
  theme_bw() 

m3 <- ggplot(as.data.frame(coord), aes(x = PC1, y = PC2, color = alpha_p)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=c(a,b)) +
  theme_bw() 

m4 <- ggplot(as.data.frame(nodeCoord), aes(x = PC1, y = PC2, color = alpha_sc)) +
  geom_point() +
  scale_colour_viridis(option='H', limits=c(a,b)) +
  theme_bw()

ggarrange(m1,m2,m3,m4, ncol=4, common.legend = T, legend="right")




#### Arc length fitting plots
load(paste0(dir, "retro_multicyclic2b_fitting.rda"))
f1 = fitting[[1]]
f2 = fitting[[2]]
f = rowMeans(do.call(rbind, fitting), na.rm=T)

pseudotime_cycle_1 = pseudotime[which(time >= 1 & time < 6)]
pseudotime_cycle_2 = pseudotime[which(time >= 6 & time < 17)]
pseudotime_lin_1 = pseudotime[which(time >= 17 & !is.na(fitting[[1]]))]
pseudotime_lin_2 = pseudotime[which(time >= 17 & !is.na(fitting[[2]]))]

# obtain arc length per region of trajectory
arclen_cycle_1 = f[which(time >= 1 & time < 6)]
arclen_cycle_2 = f[which(time >= 6 & time < 17)]
arclen_lin_1 = f1[which(time >= 17 & !is.na(fitting[[1]]))]
arclen_lin_2 = f2[which(time >= 17 & !is.na(fitting[[2]]))]

# normalizing arc lengths 
arclen_cycle_1 = arclen_cycle_1 / max(arclen_cycle_1) 
arclen_cycle_2 = (arclen_cycle_2-min(arclen_cycle_2)) / (max(arclen_cycle_2)-min(arclen_cycle_2))
arclen_lin_1 = (arclen_lin_1-min(arclen_lin_1)) / (max(arclen_lin_1)-min(arclen_lin_1))
arclen_lin_2 = (arclen_lin_2-min(arclen_lin_2)) / (max(arclen_lin_2)-min(arclen_lin_2))

# make dfs for plotting
fitting_cycle1_df = as.data.frame(cbind(pseudotime_cycle_1, arclen_cycle_1))
fitting_cycle2_df = as.data.frame(cbind(pseudotime_cycle_2, arclen_cycle_2))
fitting_lin1_df = as.data.frame(cbind(pseudotime_lin_1, arclen_lin_1))
fitting_lin2_df = as.data.frame(cbind(pseudotime_lin_2, arclen_lin_2))
colnames(fitting_cycle1_df) = colnames(fitting_cycle2_df) = colnames(fitting_lin1_df) = colnames(fitting_lin2_df) = c("pseudotime", "arc_length")

ggplot() + 
  geom_point(data = fitting_cycle1_df, aes(y=pseudotime, x=arc_length, colour="cycle-1")) + 
  geom_point(data = fitting_cycle2_df, aes(y=pseudotime, x=arc_length, colour="cycle-2")) + 
  geom_point(data = fitting_lin1_df, aes(y=pseudotime, x=arc_length, colour="lin-1")) + 
  geom_point(data = fitting_lin2_df, aes(y=pseudotime, x=arc_length, colour="lin-2")) + 
  theme_bw() + 
  labs(colour="membership") +
  xlab("normalized arc length") + 
  theme(legend.position = "bottom")



