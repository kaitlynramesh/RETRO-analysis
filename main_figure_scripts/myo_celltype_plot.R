# ??
# FIGURES (??)

# Load necessary libraries
library(ggplot2)
library(dplyr)

stack_area_plot <- function(cell_type, pseudotime, nbreaks) {
  
  # Combine the data into a dataframe
  data <- data.frame("cell_type" = cell_type[order(pseudotime)], 
                     "pseudotime" = sort(pseudotime))
  
  # Bin the pseudotime into 10 bins (you can adjust the number of bins)
  data$pseudotime_bin <- cut(data$pseudotime, breaks = nbreaks, labels = FALSE)
  
  # Calculate proportions for each cell type within each pseudotime bin
  # proportions <- data %>%
  #   group_by(pseudotime_bin, cell_type) %>%
  #   summarise(count = n()) %>%
  #   ungroup() %>%
  #   group_by(pseudotime_bin) %>%
  #   mutate(proportion = count / sum(count)) %>%
  #   ungroup()
  
  proportions <- data %>%
    dplyr::group_by(pseudotime_bin, cell_type) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(pseudotime_bin) %>%
    dplyr::mutate(proportion = count / sum(count)) %>%
    dplyr::ungroup()
  
  # Add the pseudotime values back to the proportions object
  proportions$pseudotime <- data$pseudotime[match(proportions$pseudotime_bin, data$pseudotime_bin)]
  
  # Create the stacked area plot
  ggplot(proportions, aes(x = pseudotime, y = proportion, fill = cell_type)) +
    geom_area(stat = "identity", position = "stack") +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal() +
    labs(x = "Pseudotime", y = "Proportion", fill = "Cell Type") +
    theme(legend.position = "right") + 
    scale_x_continuous(limits=c(min(pseudotime),max(pseudotime)))  # Adjust x-axis ticks for pseudotime
}

# load datasets
myo_url = "https://raw.githubusercontent.com/kaitlynramesh/RETRO-analysis/main/real/scd_myo.rds"
scd_myo = readRDS(gzcon(url(myo_url)))

# load RETRO pseudotime object
dir = "~/KaitlynRRStudio/RETRO-analysis/benchmark/"
retro_pt_obj = readRDS(file=paste0(dir, "RETRO_PT_myo.rds"))

cell_type = scd_myo@phenoData@data[["cell_type"]]
pseudotime = retro_pt_obj@pseudotime

stack_area_plot(cell_type, pseudotime, nbreaks=5)
