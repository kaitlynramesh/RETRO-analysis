# Infer kinetics distinguishing between developmental trajectories in myelopoiesis   
# Figures (??)

library(viridis)
library(ggplot2)
library(pracma)
library(ggpubr)

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

# Load data
myelo_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-data/main/real/scd_myelo.rds"
scd_myelo = readRDS(gzcon(url(myelo_url)))

# load pseudotime/fitting object


# parameters for analysis
coord = scd_myelo7$UMAP # umap projection
nl = 4 # num lineages
r = 2
kinetic_mat = matrix(nrow = nrow(coord), ncol = nl)

for(l in seq(nl)) {
  coord_subset = coord[lin_membership[[l]],]
  pt_subset = pseudotime[lin_membership[[l]]]
  
  alpha = calculate_alpha(coord_subset, pt_subset, time, radius=r)
  
  kinetic_mat[lin_membership[[l]],l] = alpha 
}

retro_kinetic = rowMeans(kinetic_mat, na.rm=T)

cell_type = as.factor(as.character(scd_myelo7$Cell_Type))
abbrev_celltype = c("AdM", "EmM", "HEn", "iPSC", "Mac", "Mast", "Mega", "MDP", "Mon", "MP", "NMP", "PS")
levels(cell_type) <- abbrev_celltype # levels(cell_type)[order(abbrev_celltype)]

colnames(coord) = c("UMAP1", "UMAP2")
g_ct = ggplot(data=as.data.frame(coord), aes(x=UMAP1, y=UMAP2, colour=cell_type)) + 
  geom_point() + 
  theme_bw()
g_k = ggplot(data=as.data.frame(coord), aes(x=UMAP1, y=UMAP2, colour=retro_kinetic)) + 
  scale_colour_viridis(option='H') + 
  geom_point() + 
  theme_bw() + theme(legend.position = "right")

ggarrange(g_ct, g_k)

psupertime_density(umap, scd_myelo7$Time, pseudotime)

ggplot(data=as.data.frame(umap), aes(x=UMAP1, y=UMAP2, colour=pseudotime)) + 
  geom_point() + 
  scale_colour_continuous(low="blue", high="orange") + 
  theme_bw()



