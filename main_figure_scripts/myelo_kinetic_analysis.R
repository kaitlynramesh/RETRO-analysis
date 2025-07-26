# Investigate kinetics of developmental trajectories in myelopoiesis based on inferred pseudotime  
# Figures (??) and (??)

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


#### Myelopoiesis data

# load data
myelo_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/real/scd_myelo.rds"
scd_myelo = readRDS(gzcon(url(myelo_url)))

# load pseudotime/fitting object
dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/" ## fix this!!
myelo_pt_obj = readRDS(file=paste0(dir, "RETRO_PT_myelo.rds"))

time = myelo_pt_obj@time
cell_type = scd_myelo@phenoData@data[["cell_type"]]
coord = scd_myelo@experimentData@other[["UMAP"]]
lin_membership = myelo_pt_obj@lin_membership
pseudotime = myelo_pt_obj@pseudotime

nl = 4 # num lineages
r = 2 # radius to determine neighbors
kinetic_mat = matrix(nrow = nrow(coord), ncol = nl) 

for(l in seq(nl)) { 
  coord_subset = coord[lin_membership[[l]],]
  pt_subset = pseudotime[lin_membership[[l]]]
  
  alpha = calculate_alpha(coord_subset, pt_subset, time, radius=r)
  
  kinetic_mat[lin_membership[[l]],l] = alpha 
}
retro_kinetic = rowMeans(kinetic_mat, na.rm=T) # average kinetic value across lineages

# plotting (1) cell type and (2) kinetic changes
colnames(coord) = c("UMAP1", "UMAP2")
ggplot(data=as.data.frame(coord), aes(x=UMAP1, y=UMAP2, colour=cell_type)) + 
  geom_point() + 
  theme_bw()

ggplot(data=as.data.frame(coord), aes(x=UMAP1, y=UMAP2, colour=retro_kinetic)) + 
  scale_colour_viridis(option='H') + 
  geom_point() + 
  theme_bw() + theme(legend.position = "right")

ggplot(data=as.data.frame(coord), aes(x=UMAP1, y=UMAP2, colour=cell_type)) + 
  geom_point() + 
  theme_bw()
