## compute occupancy across Eastern Cordillera

## housekeeping ----
compute_posterior = TRUE
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
    dt_ranges <- lapply(fnames, function(x) data.table(readRDS(x))) %>%
        setNames(species_names) %>%
        rbindlist(., idcol = "species") %>%
        setNames(c("species", "id_cell"))
    
    # indexing into dt_cov ultimately saves 1 numeric column (so 2 in cols 
    # equivalent) after ele gets removed
    setkey(dt_ranges, id_cell)
    setkey(dt_cov, id_cell)
    dt_cov[, row_id := 1:.N]
    
    dt_ranges[dt_cov, `:=`(cov_index = i.row_id)]
    
    # There are a much smaller number of unique elevations than unique cell ids,
    # and can therefore be much more memory efficient by mapping from this 
    # unique set of elevations to relevs, and then indexing these into the main
    # DT
    
    # create DT with all combinations of species and unique elevations
    dt_relev <- CJ(species = species_subset,
                   ele = unique(dt_cov[,ele]))
    
    # calculate relevs in this DT
    setkey(dt_relev, species)
    setkey(dt_ele_scaling, species)
    dt_relev[dt_ele_scaling, 
                  `:=`(relev = (ele - elev_median)/(elev_breadth/2),
                       id_row = 1:.N)]
    
    # in order to create the index into this relev DT, need to briefly add an
    # ele column to dt_ranges, create the index, and then can remove. 
    # note: the advantage of having a column that indexes into dt_ele_scaled, 
    # rather than simply having a relev column, is that the index is int, while
    # the latter is numeric (i.e. net saving of 1 int value column)
    dt_ranges[,ele := dt_cov$ele[cov_index]]
    
    setkey(dt_relev, species, ele)
    setkey(dt_ranges, species, ele)
    dt_ranges[dt_relev, relev_index := i.id_row]
    dt_ranges[, ele := NULL]
    
    ## lastly, want to create the index to index into dt_coef
    setkey(dt_ranges, species)
    setkey(dt_coef, species)
    dt_coef[,row_index := 1:.N]
    dt_ranges[dt_coef[draw == 1], coef_index := i.row_index]
    
    # finally, create cluster lookup DT. 
    # note: indexing into this has a net saving of 1 column
    cluster_lookup <- unique(dt_ranges[,.(id_cell)])
    cluster_lookup[,`:=`(id_cl = clump_id_cell(id_cell, 2),
                         id_sr = clump_id_cell(id_cell, 100),
                         row_id = 1:.N)]
    
    
    setkey(cluster_lookup, id_cell)
    setkey(dt_ranges, id_cell)
    dt_ranges[cluster_lookup, cluster_index := i.row_id]
    dt_seff <- unique(dt_ranges[cluster_lookup,.(species, id_cl, id_sr)])
    dt_seff[, eff_cl_sp := rnorm(.N, 0, dt_spatial[1, sd_sp_cl]), keyby = .(species)]
    dt_seff[, eff_sr_sp := rnorm(.N, 0, dt_spatial[1, sd_sp_cl]), keyby = .(species, id_sr)]
    
    max_cl <- max(cluster_lookup$id_cl)
    max_sr <- max(cluster_lookup$id_sr)
    
    ## calculate occupancy ----
    print("starting occupancy calculation:")
    print(Sys.time() - start)
    
   
    peakRAM::peakRAM({
        start <- Sys.time()
    dt_ranges[,occ := dt_coef$lpo[coef_index]]
    dt_ranges[, occ := occ + 
                  dt_coef$relev_term[coef_index] * dt_relev$relev[relev_index]]
    
    dt_ranges[, occ := occ +
                  dt_coef$relev2_term[coef_index] * dt_relev$relev[relev_index]^2]
    
    dt_ranges[, occ := occ +
                  rnorm(max_sr, 0, dt_spatial[1, sd_sr])[cluster_lookup$id_sr[cluster_index]]]
    
    dt_ranges[,`:=`(
        occ = occ + 
            # pasture
            (1 - dt_cov$fc[cov_index]/100) * (dt_coef$pasture[coef_index] * 1))] 
    
    dt_ranges[, occ := occ + 
            (1 - dt_cov$fc[cov_index]/100) * (dt_coef$woody_veg_sc_term[coef_index] * dt_cov$woody_veg_sc[cov_index]) 
            ]
    
    dt_ranges[, occ := occ + # forest
            (dt_cov$fc[cov_index]/100) * (dt_coef$pasture[coef_index] * -1)]
    
    dt_ranges[, occ := occ + 
            (dt_cov$fc[cov_index]/100) * (dt_coef$woody_veg_sc_term[coef_index] * (1 - woody_veg_cent)/woody_veg_scale)
        ]
    
    # rnorm(max_cl, 0, dt_spatial[i, sd_sp_cl])[cluster_lookup$id_cl[dt_ranges$cluster_index]]
    setkey(dt_ranges, id_cell)
    setkey(cluster_lookup, id_cell)
    
    # dt_ranges[, occ := boot::inv.logit(
    #     occ +
    #         dt_seff[cluster_lookup[dt_ranges, id_cl]] +
    #         rnorm(max_sr, 0, dt_spatial[1, sd_sp_sr])[cluster_lookup[dt_ranges, id_sr]]
    #     )]
    setkey(dt_seff, species, id_cl)
    # setkey(dt_ranges, species, cluster_lookup$id_cl[cluster_index])
    # dt_ranges[, occ := 
    #     occ +
    #         dt_seff[cluster_lookup[dt_ranges, id_cl]] +
    #         rnorm(max_sr, 0, dt_spatial[1, sd_sp_sr])[cluster_lookup[dt_ranges, id_sr]]
    # )]
    
    print(Sys.time() - start)
    })

    # need to work out the numbers of clusters by species
    dt_ranges$id_cell
    
    setkey(dt_ranges, id_cell)
    setkey(cluster_lookup, id_cell)
    cluster_lookup[dt_ranges][,sum, by=species]
    cluster_lookup
    
    object.size(dt_ranges)/1e6
    rep(c(1, 2, 3), c(5, 5, 1))
    
    for(i in 1:28) {
        print(i)
        peakRAM::peakRAM({
            start <- Sys.time()
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
                      rnorm(max_sr, 0, dt_spatial[i, sd_sr])[id_sr]]
        
        # note rnorm may be over-allocated for some species. 
        dt_ranges[, p := boot::inv.logit(occ +
                                             rnorm(max_cl, 0, dt_spatial[i, sd_sp_cl])[id_cl] + 
                                             rnorm(max_sr, 0, dt_spatial[i, sd_sp_sr])[id_sr]), 
                  keyby = "species"]
        
        saveRDS(dt_ranges$p, paste0("outputs/posterior_", i, ".rds"))
        print(Sys.time() - start)
        })
        
        Sys.time() - start;
        object.size(dt_ranges)/1e6
        1595
        #set key: 1592/1598
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
    
    
    lc <- readRDS("outputs/input_layers_stars.rds")
    i <- "Zonotrichia_capensis" #species_list[1]
    lc$species_i <- ifelse(is.na(lc$ele), NA, 0)
    lc$species_i[dt_ranges[species == i, id_cell]] <-
        dt_ranges[species == i, occ]
    
    plot(lc["species_i"], breaks = "equal", 
         col=viridis::inferno(20), downsample = 4,
         nbreaks=21)
    library(zoom)
    zm()
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
