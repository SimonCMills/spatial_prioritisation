## compute occupancy across Eastern Cordillera

## housekeeping ----
library(data.table); library(dplyr);

# function to clump cell ids (for including spatial terms)
clump_id_cell <- function(id_cell, grain) {
    x_pos <- ceiling(id_cell/3031)
    y_pos <- id_cell %% 3031 
    
    y_new <- ceiling(y_pos/grain)
    x_new <- ceiling(x_pos/grain)
    (y_new - 1) * max(x_new) + x_new
}

# extract some scaling parameters for use downstream 
# note: this file was originally saved at: 
# stan_out/EC_model_run6_woody_veg_effs/cu_cov.rds
cu_cov <- readRDS("outputs/cu_cov.rds")
woody_veg_cent <- attr(cu_cov$woody_veg_sc, "scaled:center")
woody_veg_scale <- attr(cu_cov$woody_veg_sc, "scaled:scale")
dt_ele_scaling <- cu_cov[, c("species", "elev_breadth", "elev_median")] %>%
    unique %>%
    as.data.table
rm(cu_cov)


start <- Sys.time()

## read in range info ----
fnames <- list.files("outputs/in_range/", full.names = TRUE)
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
coefs <- readRDS("outputs/lpo_and_coefs_EC.rds")
# for a first glance output, summarise across rows (30 iterations). Obviously
# will later change this to keep full posterior through analysis pipeline. 
dt_coef <- data.table(species = rep(coefs$species, 28), # recycled
                      draw = rep(1:28, each = length(coefs$species)),
                      relev_term = as.vector(coefs$relev_term), 
                      relev2_term = as.vector(coefs$relev2_term), 
                      woody_veg_sc_term = as.vector(coefs$woody_veg_sc_term),
                      pasture = as.vector(coefs$pasture),
                      lpo = as.vector(coefs$lpo))

dt_spatial <- coefs$spatial_terms


## join range with cov DTs & scale elev ----
dt_ranges[dt_cov, `:=`(ele = i.ele, 
                       fc=i.fc, 
                       woody_veg_sc = i.woody_veg_sc), on="id_cell"]
dt_ranges[dt_ele_scaling, 
          relev := (ele - elev_median)/(elev_breadth/2), 
          on="species"]

dt_ranges[,ele := NULL]

## calculate occupancy ----
print("starting occupancy calculation:")
print(Sys.time() - start)

## 20 km2 is 113 * 113 cells
dt_ranges[,`:=`(occ = 0,
                p = 0,
                id_cl = clump_id_cell(id_cell, 2),
                id_sr = clump_id_cell(id_cell, 100))]

max_cl <- max(dt_ranges$id_cl)
max_sr <- max(dt_ranges$id_sr)


for(i in 1:28) {
    print(i)
    print(paste0("elapsed time: ", Sys.time() - start))
    
    dt_ranges[dt_coef[draw == i,], 
              occ := i.lpo + 
                  i.relev_term * relev +
                  i.relev2_term * relev^2 +
                  # pasture
                  (1 - fc/100) * (i.pasture * 1) + 
                  (1 - fc/100) * (i.woody_veg_sc_term * woody_veg_sc) + 
                  # forest
                  (fc/100) * (i.pasture * -1) +
                  (fc/100) * (i.woody_veg_sc_term * (1 - woody_veg_cent)/woody_veg_scale) + 
                  # spatial terms
                  rnorm(max_sr, 0, dt_spatial[i, sd_sr])[id_sr], 
              on="species"]
    
    # note rnorm may be over-allocated for some species. 
    dt_ranges[, p := boot::inv.logit(occ +
        rnorm(max_cl, 0, dt_spatial[i, sd_sp_cl])[id_cl] + 
        rnorm(max_sr, 0, dt_spatial[i, sd_sp_sr])[id_sr]), 
        by = "species"]
    
    saveRDS(dt_ranges$p, paste0("outputs/posterior_", i, ".rds"))
}
