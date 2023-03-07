## compute occupancy across Eastern Cordillera

## housekeeping ----
compute_posterior = FALSE
compute_summary = FALSE
plot_figures = TRUE

library(data.table); library(dplyr);

# function to clump cell ids (for including spatial terms)
clump_id_cell <- function(id_cell, grain) {
    x_pos <- ceiling(id_cell/3031)
    y_pos <- id_cell %% 3031 
    
    y_new <- ceiling(y_pos/grain)
    x_new <- ceiling(x_pos/grain)
    (y_new - 1) * max(x_new) + x_new
}

if(compute_posterior) {
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
    species_subset <- c("Zonotrichia_capensis", 
                        "Sporophila_luctuosa",
                        "Aburria_aburri",
                        "Sericossypha_albocristata",
                        "Pyrrhomyias_cinnamomeus",
                        "Odontophorus_strophium")
    fnames <- fnames[grepl(paste0(species_subset, collapse="|"), fnames)]
    
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
    saveRDS(dt_ranges[,list(species, id_cell, ele, relev, fc)], 
            "outputs/dt_ranges_example_species.rds")
    dt_ranges[,ele := NULL]
    
    ## calculate occupancy ----
    print("starting occupancy calculation:")
    print(Sys.time() - start)
    
    ## 20 km2 is 113 * 113 cells
    # okay, will actually reduce it by a resolution of 4 (i.e. 16 times smaller)
    dt_ranges[,`:=`(occ = 0,
                    p = 0,
                    id_cl = clump_id_cell(id_cell, 2),
                    id_sr = clump_id_cell(id_cell, 100),
                    id_store = clump_id_cell(id_cell, 4))]
    
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
}    

if(compute_summary) {
    rm("dt_ranges", "dt_cov")
    fnames <- list.files("outputs/", "posterior_", full.names=TRUE)
    dat <- lapply(fnames, readRDS)
    df <- as.matrix(bind_cols(dat))
    rm(dat)
    
    summ <- data.table(estimate = matrixStats::rowMeans2(df), 
                       logit_sd = matrixStats::rowSds(boot::logit(df)))
    saveRDS(summ, "outputs/5species_summary.rds")
}

if(plot_figures) {
    
}
# lwr = matrixStats::rowQuantiles(df, probs = .25),
#                       upr = matrixStats::rowQuantiles(df, probs = .75))


# what are the quantities that we want?
# (1) the distribution of sum occupancy 

if(plot_figures) {
    library(ggplot2); library(stars)
    lc <- readRDS("outputs/input_layers_stars.rds")
    summ <- readRDS("outputs/5species_summary.rds")
    dt_ranges <- readRDS("outputs/dt_ranges_example_species.rds")
    species_list <- unique(dt_ranges$species)
    
    dt_ranges[,estimate := summ$estimate]
    dt_ranges[,sd := summ$log_sd]
    rm(summ)
    
    for(i in species_list) {
        dt_ranges[species == i] %>%
            .[sample(1:.N, 1000)] %>%
            ggplot(aes(ele, estimate, col = fc)) + 
            geom_point() +
            scale_colour_viridis_c() +
            theme_classic() +
            theme(panel.grid = element_blank(),
                  panel.background = element_rect(colour="black"),
                  axis.text = element_text(colour = "black")) +
            labs(x = "Elevation", 
                 y = "Pr(occupancy)",
                 colour = "Forest\ncover") 
        ggsave(paste0("figures/occupancy_elevation_", i, ".png"), 
               width = 70 * 2, 
               height = 40 * 2, 
               dpi = 400,
               units="mm")
    }
    
    i <- species_list[1]
    lc$species_i <- ifelse(is.na(lc$ele), NA, 0)
    lc$species_i[dt_ranges[species == i, id_cell]] <-
        dt_ranges[species == i, estimate]
    
    plot(lc["species_i"], breaks = "equal", col=viridis::cividis(20), nbreaks=21)
    
    
    # get Colombia bounding box (this will be the extent of the inset plot)
    CO_buffer <- rnaturalearth::ne_countries(country="Colombia") %>%
        st_as_sf %>%
        st_transform(., crs=flat) %>%
        st_buffer(50000) %>%
        st_bbox()
    
    # get country outlines
    country_borders <- rnaturalearth::ne_countries(scale=10) %>%
        st_as_sf %>%
        st_crop(., bbox_inset) %>%
        st_transform(., crs=flat) %>% 
        st_crop(., CO_buffer)
    
    # get mountain ranges
    mountains <- readRDS("../range_shifts/data/mountain_polygons.rds") %>% 
        st_crop(., lc)
    
    # plot
    for(i in species_list) {
        lc$species_i <- ifelse(is.na(lc$ele), NA, 0)
        lc$species_i[dt_ranges[species == i, id_cell]] <-
            dt_ranges[species == i, estimate]
        
        png(paste0("figures/map_occ_", i, ".png"), res = 300, height = 200, 
            width = 150, units="mm")
        plot(lc["species_i"], 
             breaks = "equal", 
             col=viridis::inferno(50), 
             nbreaks=51, 
             downsample = 20,
             reset = FALSE,
             main = "")
        
        plot(mountains[,1], col="grey90", border="grey90", main="", add = T)
        
        plot(lc["species_i"], 
             breaks = "equal", 
             col=viridis::inferno(50), 
             nbreaks=51, 
             downsample = 3,
             reset = FALSE,
             main = "",
             add = T)
        plot(st_as_sfc(st_bbox(lc)), add=T)
        dev.off()
    }
}
