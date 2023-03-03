# convert range polygons to data.table of cell ids where each species is present
# Now doing this locally, because can't face faffing with spatial libraries. 
# Transfer these files across to HPC after and HPC from there on out
# Notes:
# - ranges aren't buffered.. should they be?
# - run on HPC
# packages
library(sf); library(stars); library(data.table); library(dplyr)

# script for cluster

# 821 species have a match. Were are the others?
# TO DO: resolve taxonomic mismatches
sum(ayerbe_maps$Species %in% gsub("_", " ", species_list))

ayerbe_maps_trim <- ayerbe_maps %>%
    filter(Species %in% species_list)

species_subset <- gsub(" ", "_", ayerbe_maps_trim$Species)

## sod it just do it locally. 
species_subset <- gsub(" ", "_", ayerbe_maps_trim$Species)

for(i in 167:length(species_subset)) {
    print(i)
    if(i %% 50 == 0) file.remove(list.files(tempdir(), full.names = TRUE)) 
    
    r <- ayerbe_maps_trim[i,] %>% 
        mutate(in_range = 1) %>%
        st_rasterize(., template = grid)
    
    saveRDS(in_EC[r[[1]][in_EC] == 1], 
            paste0("outputs/in_range/", species_subset[i], ".rds"))
}
