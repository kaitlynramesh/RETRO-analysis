library(RETRO)
library(Biobase) # data storage + processing
library(doParallel)
library(parallel)

scd_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_bifurcation.rds"
scd_bifurcation = readRDS(gzcon(url(scd_url)))
print("Data downloaded.")

lambda = 1
k_range = seq(8,20,2) # set range of clusters 
num_scores = 100
filename = "RETRO_bifurcation_v51.rda"

pca = experimentData(scd_bifurcation)@other$PCA
time = phenoData(scd_bifurcation)@data[["time"]] 
time_label = as.factor(time) # factor for coloring

coordinates <- weight_coord(scd_bifurcation, weight = lambda)
coordinates_uw <- weight_coord(scd_bifurcation, weight = 0)

retro_obj = set_RETRO_class(coordinates=coordinates,
                            coordinates_uw=coordinates_uw,
                            k_range=k_range,
                            time=time)

ncores = parallel::detectCores()
registerDoParallel(cores = ncores)

retro_obj = scoring(retro_obj, k_range, num_scores)

print("RETRO done.")
save(retro_obj, file=filename)
print("File saved.")


