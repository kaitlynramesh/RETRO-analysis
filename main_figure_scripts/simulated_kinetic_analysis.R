
library(viridis)
library(ggplot2)
library(pracma)

# Function to infer alpha (kinetic value) 
calculate_alpha <- function(coord, pseudotime, time, rad_per, cyclic=FALSE, filt=FALSE) {
  
  m <- as.matrix(dist(coord)) # gene exprs. dist matrix
  radius = rad_per * max(m)
  print(radius)
  
  t = unique(time)[-c(1,length(unique(time)))] 
  breaks = sort(unique(c(min(pseudotime), t, max(pseudotime))))
  
  bins = cut(pseudotime, breaks=breaks, labels=seq(length(breaks)-1))
  bins = as.numeric(bins)
  
  neighbor_index = vector(mode="list", length=nrow(coord))
  alpha_values = seq(nrow(coord))
  delta_p_values = vector(mode="list", length=nrow(coord))
  
  for(i in 1:nrow(coord)) {
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
    neighbor_index[[i]] = index # save neighbors for smoothing step
    
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
    alpha_values[i] <- 1/norm(alpha, type='2') # save "cell-speed" 
    
    # To check linear fitting
    delta_p = pseudotime[i]-pseudotime[index] # actual ∆p
    delta_p_pred <- X %*% alpha # predicted ∆p
    delta_p_values[[i]] = list(delta_p_pred, delta_p) # save ∆p values
  }
  
  # Determined smoothed alpha values for each neighborhood
  alpha_s <- lapply(1:length(alpha_values), function(i) { 
    n = unlist(neighbor_index[[i]])
    mean_alpha <- mean(alpha_values[n])
    return(mean_alpha)
  })
  alpha_s <- log(unlist(alpha_s))

  return(list("alpha" = alpha_s, "delta_p" = delta_p_values))
}

scale_time <- function(time, pseudotime) {
  or <- as.vector(range(pseudotime)) # original range
  ar <- as.vector(range(time)) # actual range
  scaled_pt <- lapply(pseudotime, function(p) ((((p-or[1])*(ar[2]-ar[1])) / (or[2]-or[1])) + ar[1]) )
  scaled_pt <- unlist(scaled_pt)
  return(scaled_pt)
}

filter_by_cor <- function(delta_p_values, cutoff) {
  # get correlation
  cor = lapply(delta_p_values, function(x) {
    pred = as.numeric(x[[1]])
    obs = as.numeric(x[[2]])
    cor = cor.test(pred, obs)
    return(cor$estimate)})
  cor = as.numeric(cor)
  
  # get mean pred/obs ∆p values per cell
  delta_p_pred = sapply(delta_p_values, "[", 1)
  delta_p_obs = sapply(delta_p_values, "[", 2)
  delta_p_pred = sapply(delta_p_pred, mean) 
  delta_p_obs = sapply(delta_p_obs, mean)
  
  delta_p_df = cbind("delta_p_pred"=delta_p_pred,
                     "delta_p_obs"=delta_p_obs)
  # ggplot(as.data.frame(delta_p_df), aes(x=delta_p_pred,y=delta_p_obs, colour=cor)) + 
  #   labs(x="Mean observed ∆p", y="Mean predicted ∆p") + 
  #   geom_point()
  
  keep_cells = which(cor > cutoff)
  return(keep_cells)
}



##### Load data #####

bif_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_bifurcation.rds"
twocycles_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_twocycles.rds"
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
sctime = sctda_list$sctime[-length(sctda_list$sctime)] # sampling time at each cell 

r = 0.1

alpha_r <- calculate_alpha(coord, pseudotime, time, rad_per=r)
alpha_s <- calculate_alpha(coord, slingshot_pt, time, rad_per=r)
alpha_p <- calculate_alpha(coord, psupertime_pt, time, rad_per=r)
alpha_sc <- calculate_alpha(nodeCoord, scTDA_pt, sctime, rad_per=r)

cutoff = 0.5
r_cells = filter_by_cor(alpha_r$delta_p, cutoff)
s_cells = filter_by_cor(alpha_s$delta_p, cutoff)
p_cells = filter_by_cor(alpha_p$delta_p, cutoff)
sc_cells = filter_by_cor(alpha_sc$delta_p, cutoff)

alpha_r = alpha_r$alpha[r_cells]
alpha_s = alpha_s$alpha[s_cells]
alpha_p = alpha_p$alpha[p_cells]
alpha_sc = alpha_sc$alpha[sc_cells]

a <- min(alpha_r, alpha_p, alpha_s, alpha_sc)
b <- max(alpha_r, alpha_p, alpha_s, alpha_sc)
y_lim = c(a, b)

# plotting alpha
m1 <- ggplot(as.data.frame(coord[r_cells,]), aes(x = PC1, y = PC2, color = alpha_r)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") + 
  theme_bw() 

m2 <- ggplot(as.data.frame(coord[s_cells,]), aes(x = PC1, y = PC2, color = alpha_s)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") + 
  theme_bw() 

m3 <- ggplot(as.data.frame(coord[p_cells,]), aes(x = PC1, y = PC2, color = alpha_p)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") + 
  theme_bw() 

m4 <- ggplot(as.data.frame(nodeCoord[sc_cells,]), aes(x = PC1, y = PC2, color = alpha_sc)) +
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") +
  theme_bw()

ggarrange(m1,m2,m3,m4, ncol=4, common.legend = T, legend="right")


## arc-length/time fitting plot
lin_membership = retro_pt_obj@lin_membership
fitting = retro_pt_obj@arc_length

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
coord = scd_twocycles@experimentData@other[["PCA"]][["x"]] # PCA values
coord = coord[,1:2]
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

r = 0.075

alpha_r <- calculate_alpha(coord, pseudotime, time, rad_per=r, cyclic=TRUE)
alpha_s <- calculate_alpha(coord, slingshot_pt, time, rad_per=r, cyclic=TRUE)
alpha_p <- calculate_alpha(coord, psupertime_pt, time, rad_per=r, cyclic=TRUE)
alpha_sc <- calculate_alpha(nodeCoord, scTDA_pt, sctime, rad_per=r, cyclic=TRUE)

cutoff = 0.5
r_cells = filter_by_cor(alpha_r$delta_p, cutoff)
s_cells = filter_by_cor(alpha_s$delta_p, cutoff)
p_cells = filter_by_cor(alpha_p$delta_p, cutoff)
sc_cells = filter_by_cor(alpha_sc$delta_p, cutoff)

alpha_r = alpha_r$alpha[r_cells]
alpha_s = alpha_s$alpha[s_cells]
alpha_p = alpha_p$alpha[p_cells]
alpha_sc = alpha_sc$alpha[sc_cells]

a <- min(alpha_r, alpha_p, alpha_s, alpha_sc)
b <- max(alpha_r, alpha_p, alpha_s, alpha_sc)
y_lim = c(a, b)

# plotting alpha
m1 <- ggplot(as.data.frame(coord[r_cells,]), aes(x = PC1, y = PC2, color = alpha_r)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") + 
  theme_bw() 

m2 <- ggplot(as.data.frame(coord[s_cells,]), aes(x = PC1, y = PC2, color = alpha_s)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") +
  theme_bw() 

m3 <- ggplot(as.data.frame(coord[p_cells,]), aes(x = PC1, y = PC2, color = alpha_p)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") +
  theme_bw() 

m4 <- ggplot(as.data.frame(nodeCoord[sc_cells,]), aes(x = PC1, y = PC2, color = alpha_sc)) +
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") +
  theme_bw()

ggarrange(m1,m2,m3,m4, ncol=4, common.legend = T, legend="right")


## arc-length/time fitting plot
f = unlist(retro_pt_obj@arc_length)

pseudotime_cycle_1 = pseudotime[time <= 10]
pseudotime_cycle_2 = pseudotime[time > 10]

arclen_1 = f[time <= 10] / max(f[time <= 10]) 
arclen_2 = (f[time > 10]-min(f[time>10])) / max(f[time > 10]-min(f[time>10]))

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

r = 0.1

alpha_r <- calculate_alpha(coord, pseudotime, time, rad_per=r, cyclic=TRUE)
alpha_s <- calculate_alpha(coord, slingshot_pt, time, rad_per=r, cyclic=TRUE)
alpha_p <- calculate_alpha(coord, psupertime_pt, time, rad_per=r, cyclic=TRUE)
alpha_sc <- calculate_alpha(nodeCoord, scTDA_pt, sctime, rad_per=r, cyclic=TRUE)

cutoff = 0.5
r_cells = filter_by_cor(alpha_r$delta_p, cutoff)
s_cells = filter_by_cor(alpha_s$delta_p, cutoff)
p_cells = filter_by_cor(alpha_p$delta_p, cutoff)
sc_cells = filter_by_cor(alpha_sc$delta_p, cutoff)

alpha_r = alpha_r$alpha[r_cells]
alpha_s = alpha_s$alpha[s_cells]
alpha_p = alpha_p$alpha[p_cells]
alpha_sc = alpha_sc$alpha[sc_cells]

a <- min(alpha_r, alpha_p, alpha_s, alpha_sc)
b <- max(alpha_r, alpha_p, alpha_s, alpha_sc)
y_lim = c(a, b)

# plotting alpha
m1 <- ggplot(as.data.frame(coord[r_cells,]), aes(x = PC1, y = PC2, color = alpha_r)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") +
  theme_bw() 

m2 <- ggplot(as.data.frame(coord[s_cells,]), aes(x = PC1, y = PC2, color = alpha_s)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") +
  theme_bw() 

m3 <- ggplot(as.data.frame(coord[p_cells,]), aes(x = PC1, y = PC2, color = alpha_p)) + 
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") +
  theme_bw() 

m4 <- ggplot(as.data.frame(nodeCoord[sc_cells,]), aes(x = PC1, y = PC2, color = alpha_sc)) +
  geom_point() +
  scale_colour_viridis(option='H', limits=y_lim, name="alpha") +
  theme_bw()

ggarrange(m1,m2,m3,m4, ncol=4, common.legend = T, legend="right")




#### Arc length fitting plots
fitting = retro_pt_obj@arc_length
f1 = fitting[[1]]
f2 = fitting[[2]]
f = rowMeans(do.call(cbind, fitting), na.rm=T)

pseudotime_cycle_1 = pseudotime[which(time >= 1 & time < 6)]
pseudotime_cycle_2 = pseudotime[which(time >= 6 & time < 14)]
pseudotime_lin_1 = pseudotime[which(time >= 14 & !is.na(fitting[[1]]))]
pseudotime_lin_2 = pseudotime[which(time >= 14 & !is.na(fitting[[2]]))]

# obtain arc length per region of trajectory
arclen_cycle_1 = f[which(time >= 1 & time < 6)]
arclen_cycle_2 = f[which(time >= 6 & time < 14)]
arclen_lin_1 = f1[which(time >= 14 & !is.na(fitting[[1]]))]
arclen_lin_2 = f2[which(time >= 14 & !is.na(fitting[[2]]))]

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



