## cluster admin: copy files across to cluster for running there. 

## housekeeping ----
library(sf); library(stars); library(data.table); library(dplyr)

# files
ayerbe_maps <- readRDS("inputs/initial_ayerbe_map_sf.RDS")
species_list <- readRDS("stan_out/EC_model_run6_woody_veg_effs/fd.rds")$data %>%
    pull(species) %>%
    unique %>%
    gsub("_", " ", .)

grid <- read_stars("outputs/proportion_unmasked_177m.tif") %>%
    rename(grd = 1) %>%
    mutate(grd = 0)

in_EC <- readRDS("outputs/input_layers_dt.rds")$id_cell

# write to cluster
# save(file = "Z:/edwards_lab1/User/bo1scm/spatial_prioritisation/inputs/files_for_cluster.rdata",
#      list = c("ayerbe_maps",
#               "species_list",
#               "grid",
#               "in_EC"))

file.copy("code/map_range_polygons_to_cell_id.R",
          "Z:/edwards_lab1/User/bo1scm/spatial_prioritisation/code/", TRUE)
