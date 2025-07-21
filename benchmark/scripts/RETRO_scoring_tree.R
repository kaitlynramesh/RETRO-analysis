library(RETRO)
library(Biobase) # data storage + processing
library(doParallel)
library(parallel)

scd_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_tree.rds"
scd_tree = readRDS(gzcon(url(scd_url)))
print("Data downloaded.")

lambda = 0
k_range = seq(8,20,2)
num_scores = 100
threshold = 0.1
filename = "RETRO_tree_v51.rda"

pca = experimentData(scd_tree)@other$PCA
time = phenoData(scd_tree)@data[["time"]] 
time_label = as.factor(time) # factor for coloring

coordinates <- weight_coord(scd_tree, weight = lambda)
coordinates_uw <- weight_coord(scd_tree, weight = 0)

retro_obj = set_RETRO_class(coordinates=coordinates,
                            coordinates_uw=coordinates_uw,
                            k_range=k_range,
                            time=time,
                            threshold=threshold)
                            
ncores = parallel::detectCores()
registerDoParallel(cores = ncores)

retro_obj = scoring(retro_obj, k_range, num_scores)

print("RETRO done.")
save(retro_obj, file=filename)
print("File saved.")


