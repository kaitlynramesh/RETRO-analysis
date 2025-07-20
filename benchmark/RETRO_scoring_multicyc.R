library(RETRO)
library(Biobase) # data storage + processing
library(doParallel)
library(parallel)

scd_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/synthetic/scd_multicycles.rds"
scd_multicyc = readRDS(gzcon(url(scd_url)))
print("Data downloaded.")

lambda = 1
k_range = seq(8,20,2) # set range of clusters 
num_scores = 100
period = 4
filename = "RETRO_multicycles_v51.rda"

pca = experimentData(scd_multicyc)@other$PCA
time = phenoData(scd_multicyc)@data[["time"]] 
time_label = as.factor(time) # factor for coloring

coordinates <- weight_coord(scd_multicyc, weight = lambda)
coordinates_uw <- weight_coord(scd_multicyc, weight = 0)

retro_obj = set_RETRO_class(coordinates=coordinates,
                            coordinates_uw=coordinates_uw,
                            k_range=k_range,
                            time=time,
                            period=period)

ncores = parallel::detectCores()
registerDoParallel(cores = ncores)

retro_obj = scoring(retro_obj, k_range, num_scores)

print("RETRO done.")
save(retro_obj, file=filename)
print("File saved.")


