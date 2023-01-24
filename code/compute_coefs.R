# extract linear predictor for non-spatially varying component of model

# initially start with junk model (that has a similar construction)

library(flocker); library(dplyr); library(brms); library(data.table)

fit_full <- readRDS("../sparing_sharing_elevation/warmup_27-10-2022_1032_j2789255.rds")
fit <- fit_full
thinning <- 5
for(i in seq_len(3)) {
    fit$fit@sim$samples[[i]] <- fit_full$fit@sim$samples[[i]][seq(1, 500, thinning),]
    attr(fit$fit@sim$samples[[i]], "sampler_params") <-
        attr(fit_full$fit@sim$samples[[i]], "sampler_params")[seq(1, 500, thinning),]
}

fit$fit@sim$thin <- thinning
fit$fit@sim$n_save <- rep(100, 3)
rm(fit_full)

# get closure unit covariates 
cu_cov <- readRDS("../sparing_sharing_elevation/cu_cov.rds")

uc_species <- cu_cov %>%
    select(species, 
           phylo, 
           elev_median, 
           elev_median_sc, 
           elev_breadth, 
           elev_breadth_sc, 
           forest_dep, 
           mass, 
           log_mass_sc, 
           diet_inv, 
           diet_carn, 
           diet_fruitnect, 
           diet_gran, 
           strata_index, 
           strata_aerial) %>%
    unique()


uc_full <- uc_species %>%
    mutate(not_forest = 0, 
           woody_veg_sc = 0, 
           relev = 0, 
           relev2 = 0, 
           observer_pt = "SCM", 
           paramo = -1) %>%
    filter(!duplicated(species)) # drop the second Grallaria_rufula

uc_full2 <- bind_rows(uc_full %>% mutate(not_forest = -1), 
                      uc_full %>% mutate(not_forest = 1))

ec <- list(
    time = matrix(0, nrow(uc_full2), 4)
)

obs <- matrix(0, nrow(uc_full2), 4)

newdata <- flocker::make_flocker_data(obs, uc_full2, ec)

lpo <- t(
    brms::posterior_linpred(
        fit, 
        dpar = "occ",
        ##### TODO: make sure that when fm3 is multi-chain, the draw_ids is getting out the same
        ##### iterations as iter does elsewhere!
        
        newdata = newdata$data[1:nrow(uc_full2), ], 
        re_formula = ~ (1 + not_forest + paramo + woody_veg_sc + relev + relev2 | species) +
            (1 + not_forest | gr(phylo, cov = A)
            )
    )
)

uc_index <- 1:nrow(uc_full)
template_mat <- matrix(1, nrow = length(uc_index), ncol = ncol(lpo))

## now get out relev and woody veg coefs
relev_term1 <- template_mat %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_relev", permuted = FALSE))) 

relev_term2_raw <- rstan::extract(fit$fit, 
                                  paste0("r_species__occ[", uc_full2$species[uc_index], ",relev]"), 
                                  permuted = FALSE)

relev_term2 <- cbind(t(relev_term2_raw[,1,]), 
                     t(relev_term2_raw[,2,]), 
                     t(relev_term2_raw[,3,]))

relev_term <- relev_term1 + relev_term2 


## relev2 
relev2_term1 <- template_mat %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_relev2", permuted = FALSE))) 

relev2_term2_raw <- rstan::extract(fit$fit, 
                                  paste0("r_species__occ[", uc_full2$species[uc_index], ",relev2]"), 
                                  permuted = FALSE)

relev2_term2 <- cbind(t(relev2_term2_raw[,1,]), 
                      t(relev2_term2_raw[,2,]), 
                      t(relev2_term2_raw[,3,]))

relev2_term <- relev2_term1 + relev2_term2 

## woody_veg
unique(uc_full$forest_dep)

template_mat_HD <- matrix(uc_full$forest_dep == "High", 
                          nrow = length(uc_index), ncol = ncol(lpo))
template_mat_MD <- matrix(uc_full$forest_dep == "Medium", 
                          nrow = length(uc_index), ncol = ncol(lpo))
template_mat_LD <- matrix(uc_full$forest_dep == "Low", 
                          nrow = length(uc_index), ncol = ncol(lpo))

woody_veg_sc_term1 <- template_mat %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_woody_veg_sc", permuted = FALSE))) +
    template_mat_HD %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_forest_depHigh:woody_veg_sc", permuted = FALSE))) +
    template_mat_MD %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_forest_depMedium:woody_veg_sc", permuted = FALSE))) +
    template_mat_LD %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_forest_depLow:woody_veg_sc", permuted = FALSE)))

woody_veg_sc_term2_raw <- rstan::extract(fit$fit, 
                                         paste0("r_species__occ[", uc_full2$species[uc_index], ",woody_veg_sc]"), 
                                         permuted = FALSE)

woody_veg_sc_term2 <- cbind(t(woody_veg_sc_term2_raw[,1,]), 
                            t(woody_veg_sc_term2_raw[,2,]), 
                            t(woody_veg_sc_term2_raw[,3,]))

woody_veg_sc_term <- woody_veg_sc_term1 + woody_veg_sc_term2

# forest then pasture 
out <- list(species = uc_full$species[uc_index],
            lpo_forest = lpo[1:(nrow(lpo)/2),],
            lpo_pasture = lpo[(nrow(lpo)/2 + 1):nrow(lpo),],
            relev_term = relev_term,
            relev2_term = relev2_term,
            woody_veg_sc_term = woody_veg_sc_term)

saveRDS(out, "outputs/lpo_and_coefs_EC.rds")

## spatially varying terms ----
spatial_input_layers <- readRDS("outputs/spatial_input_layers.rds")
generic_input_layers <- readRDS("outputs/generic_input_layers.rds")

spatial_input_layers[, `:=`(elev_median = NULL, 
                            elev_breadth = NULL)]

# get scaling for woody vegetation
wv_centering <- attr(cu_cov$woody_veg_sc, "scaled:center")
wv_scaling <- attr(cu_cov$woody_veg_sc, "scaled:scale")
wv_val_forest <- cu_cov %>% 
    filter(not_forest == -1) %>%
    pull(woody_veg_sc) %>%
    unique %>% as.vector

generic_input_layers[, pct_cell_pasture := pct_unmasked * (1-pct_forest)]
generic_input_layers[, pct_cell_forest := pct_unmasked * pct_forest]
generic_input_layers[, woody_veg_sc := (avg_treecover - wv_centering)/wv_scaling]

coefs <- data.table(species = out$species, 
                    relev_term = out$relev_term[,1], 
                    relev2_term = out$relev2_term[,1], 
                    woody_veg_sc_term = out$woody_veg_sc_term[,1], 
                    lpo_forest = out$lpo_forest[,1], 
                    lpo_pasture = out$lpo_pasture[,1])

spatial_input_layers[generic_input_layers, `:=`(
    woody_veg_sc = i.woody_veg_sc, 
    pct_cell_forest = i.pct_cell_forest, 
    pct_cell_pasture = i.pct_cell_pasture
), on = "index"]

spatial_input_layers[coefs, `:=`(
    logit_occ_forest = lpo_forest + woody_veg_sc_term * wv_val_forest + 
        relev_term * relev + relev2_term * relev^2, 
    lpo_pasture = i.lpo_forest, 
    logit_occ_pasture = lpo_pasture + woody_veg_sc_term * woody_veg_sc  + 
        relev_term * relev + relev2_term * relev^2), 
    on = "species"]

any(is.na(sp_i$relev))
# check <- readRDS("outputs/stars_lookup.rds") %>%
#     mutate(forest_occ = NA, pasture_occ = NA)
# sp_i <- spatial_input_layers[species == "Zonotrichia_capensis",]
sp_i <- spatial_input_layers[species == "Adelomyia_melanogenys",]

sp_i <- spatial_input_layers %>% filter(species == "Adelomyia_melanogenys")
sp_i <- spatial_input_layers %>% filter(species == "Aburria_aburri")

check[["forest_occ"]][sp_i$index] <- sp_i$woody_veg_sc#boot::inv.logit(sp_i$logit_occ_forest)
check[["pasture_occ"]][sp_i$index] <- sp_i$lpo_pasture#boot::inv.logit(sp_i$logit_occ_pasture)
# 
# 
plot(check["forest_occ"], breaks = "equal")
plot(check["pasture_occ"], breaks = "equal")

spatial_input_layers[, `:=`(
    units_forest = (pct_cell_forest * 1000^2)/(100^2*pi), 
    units_pasture = (pct_cell_pasture * 1000^2)/(100^2*pi)
    )]

spatial_input_layers[, `:=`(
    amt_occ = units_forest * boot::inv.logit(logit_occ_forest) + 
        units_pasture * boot::inv.logit(logit_occ_pasture)
)]

spatial_input_layers[,prop_range := amt_occ/sum(amt_occ), by="species"]

new <- spatial_input_layers[, sum(prop_range), by="index"]
new <- species_richness <- spatial_input_layers[, sum(amt_occ >= 1), by="index"]


check <- readRDS("outputs/stars_lookup.rds") %>%
    mutate(sum_range = NA)

check[["sum_range"]][new$index] <- new$V1
check <- check %>%
    mutate(sum_range_resc = sum_range/mean(sum_range, na.rm=T))
check
plot(check["sum_range"])
plot(check["sum_range_resc"], na.rm=T)

hist(log(check[["sum_range"]]/mean(check[["sum_range"]], na.rm=T)))
exp(2)


sp_i <- spatial_input_layers %>% 
    filter(species == "Aburria_aburri")

