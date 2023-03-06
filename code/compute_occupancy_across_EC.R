
## housekeeping ----
library(data.table); library(dplyr)

# extract some scaling parameters for use downstream 
cu_cov <- readRDS("stan_out/EC_model_run6_woody_veg_effs/cu_cov.rds")
woody_veg_cent <- attr(cu_cov$woody_veg_sc, "scaled:center")
woody_veg_scale <- attr(cu_cov$woody_veg_sc, "scaled:scale")
dt_ele_scaling <- cu_cov[, c("species", "elev_breadth", "elev_median")] %>%
    unique %>%
    as.data.table
rm(cu_cov)


start <- Sys.time()

## read in range info ----
# locally just use first 5 species (remove for full job on cluster)
fnames <- list.files("outputs/in_range/", full.names = TRUE)[1:5]
species_names <- gsub(".*in_range/(.*).rds", "\\1", fnames)

# read in
dt_ranges <- lapply(fnames, function(x) data.table(readRDS(x))) %>%
    setNames(species_names) %>%
    rbindlist(., idcol = "species") %>%
    setNames(c("species", "id_cell"))

print("reading & binding ranges done:")
print(Sys.time() - start)

## read in covariates ----
dt_cov <- readRDS("outputs/input_layers_dt.rds")

# scale woody veg values 
dt_cov[, tc := ifelse(tc > 50, 50, tc)]
# set 100% forest cover to have 0 to prevent NAs propagating downstream
dt_cov[, tc := ifelse(is.na(tc), 0, tc)]
dt_cov[, tc := (tc/100 - woody_veg_cent)/woody_veg_scale]
setnames(dt_cov, "tc", "woody_veg_sc")


## read in coefficients ----
dt_coef <- readRDS("outputs/lpo_and_coefs_EC.rds")
# for a first glance output, summarise across rows (30 iterations). Obviously
# will later change this to keep full posterior through analysis pipeline. 
dt_coef <- data.table(species = dt_coef$species, 
                      relev_term = rowMeans(dt_coef$relev_term), 
                      relev2_term = rowMeans(dt_coef$relev2_term), 
                      woody_veg_sc_term = rowMeans(dt_coef$woody_veg_sc_term),
                      pasture = rowMeans(dt_coef$pasture),
                      lpo = rowMeans(dt_coef$lpo))


## join range with cov DTs & scale elev ----
dt_ranges[dt_cov, `:=`(ele = i.ele, 
                       fc=i.fc, 
                       woody_veg_sc = i.woody_veg_sc), on="id_cell"]
dt_ranges[dt_ele_scaling, 
          relev := (ele - elev_median)/(elev_breadth/2), 
          on="species"]

dt_ranges[,ele := NULL]

## calculate occupancy ----
dt_ranges[dt_coef, occ := (
    i.lpo + 
        i.relev_term * relev +
        i.relev2_term * relev^2 +
        # pasture
        (1 - fc/100) * (i.pasture * 1) + 
        (1 - fc/100) * (i.woody_veg_sc_term * woody_veg_sc) + 
        # forest
        (fc/100) * (i.pasture * -1) +
        (fc/100) * (i.woody_veg_sc_term * (1 - woody_veg_cent)/woody_veg_scale)), 
    on="species"]


dt_ranges

## visual check 
library(stars)
readRDS("outputs/")
#[,dt_ranges[dt_coef, on="species"][
dt_ranges[dt_coef, occ := lpo_forest + i.relev_term * ele + i.relev2_term * ele^2, 
          on = "species"]

dt_ranges[species == "Aburria_aburri"]
lc$species <- NA
lc[["species"]][dt_ranges[species == "Aburria_aburri", id_cell]] <- 
    dt_ranges[species == "Aburria_aburri", occ]



function(elev, elev_median, elev_breadth)
(elev - elev_median)/elev_breadth
# what is the calculation?
occ_pasture 

dt_cov <- lc_dt