# extract linear predictor for non-spatially varying component of model
# The idea here is that there is a species-level component of the prediction 
# that we can compute the linear predictor for and only need to do once per
# (species, iteration) combination rather than repeating it for every 
# (species, cell) combination. 

# housekeeping
library(flocker); library(dplyr); library(brms); library(data.table)

fit_full <- readRDS("stan_out/EC_model_run6_woody_veg_effs/rebuilt_fit.rds")
fit <- fit_full

# thin model
# when you compute the lpo for a subset of draws it still seems to spike memory, 
# so thin down to a small object that you can work with temporarily on laptop
draws_per_chain <- ndraws(fit_full)/nchains(fit_full)
thinning <- 20
draw_ids <- seq(1, draws_per_chain, thinning)

for(i in seq_len(nchains(fit_full))) {
    fit$fit@sim$samples[[i]] <- fit_full$fit@sim$samples[[i]][draw_ids,]
    attr(fit$fit@sim$samples[[i]], "sampler_params") <-
        attr(fit_full$fit@sim$samples[[i]], "sampler_params")[draw_ids,]
}

fit$fit@sim$thin <- thinning
fit$fit@sim$n_save <- rep(length(draw_ids), nchains(fit))
rm(fit_full)

# get closure unit covariates 
# cu_cov <- readRDS("../sparing_sharing_elevation/cu_cov.rds")
fd <- readRDS("stan_out/EC_model_run6_woody_veg_effs/fd.rds")
uc_species <- fd$data %>%
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

# create flocker data object with all necessary info (note: some of this isn't
# used in generating the prediction, but still needs to be present in the object
# passed to posterior_linpred)
uc_full <- uc_species %>%
    mutate(pasture = 0, 
           woody_veg_sc = 0, 
           relev = 0, 
           relev2 = 0, 
           observer_pt = "SCM") %>%
    filter(!duplicated(species)) # drop the second Grallaria_rufula

ec <- list(
    time = matrix(0, nrow(uc_full), 4)
)

obs <- matrix(0, nrow(uc_full), 4)
newdata <- flocker::make_flocker_data(obs, uc_full, ec)

lpo <- t(
    brms::posterior_linpred(
        fit,
        dpar = "occ",
        ##### TODO: make sure that when fm3 is multi-chain, the draw_ids is getting out the same
        ##### iterations as iter does elsewhere!
        
        newdata = newdata$data[1:nrow(uc_full), ], 
        re_formula = ~ (1 + not_forest + woody_veg_sc + relev + relev2 | species) +
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
                                  paste0("r_species__occ[", uc_full$species[uc_index], ",relev]"), 
                                  permuted = FALSE)

relev_term2 <- cbind(t(relev_term2_raw[,1,]), 
                     t(relev_term2_raw[,2,]), 
                     t(relev_term2_raw[,3,]), 
                     t(relev_term2_raw[,4,]))

relev_term <- relev_term1 + relev_term2 

## relev2 
relev2_term1 <- template_mat %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_relev2", permuted = FALSE))) 

relev2_term2_raw <- rstan::extract(fit$fit, 
                                  paste0("r_species__occ[", uc_full$species[uc_index], ",relev2]"), 
                                  permuted = FALSE)

relev2_term2 <- cbind(t(relev2_term2_raw[,1,]), 
                      t(relev2_term2_raw[,2,]), 
                      t(relev2_term2_raw[,3,]),
                      t(relev2_term2_raw[,4,]))

relev2_term <- relev2_term1 + relev2_term2 

## woody_veg
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
                                         paste0("r_species__occ[", uc_full$species[uc_index], ",woody_veg_sc]"), 
                                         permuted = FALSE)

woody_veg_sc_term2 <- cbind(t(woody_veg_sc_term2_raw[,1,]), 
                            t(woody_veg_sc_term2_raw[,2,]), 
                            t(woody_veg_sc_term2_raw[,3,]),
                            t(woody_veg_sc_term2_raw[,4,]))

woody_veg_sc_term <- woody_veg_sc_term1 + woody_veg_sc_term2

## forest ----
# fixed effect terms (interactions with forest dependency)
pasture_term1 <- template_mat %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_pasture", permuted = FALSE))) +
    template_mat_HD %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_forest_depHigh:pasture", 
                                  permuted = FALSE))) +
    template_mat_MD %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_forest_depMedium:pasture", 
                                  permuted = FALSE))) +
    template_mat_LD %*%
    diag(as.vector(rstan::extract(fit$fit, "b_occ_forest_depLow:pasture", 
                                  permuted = FALSE)))

# independent random effects
pasture_term2_raw <- rstan::extract(fit$fit, 
                                    paste0("r_species__occ[", uc_full$species[uc_index], ",pasture]"), 
                                    permuted = FALSE)

pasture_term2 <- cbind(t(pasture_term2_raw[,1,]), 
                       t(pasture_term2_raw[,2,]), 
                       t(pasture_term2_raw[,3,]),
                       t(pasture_term2_raw[,4,]))

# phylogenetic effects
pasture_term3_raw <- rstan::extract(fit$fit, 
                                    paste0("r_phylo__occ[", uc_full$phylo[uc_index], ",pasture]"), 
                                    permuted = FALSE)

pasture_term3 <- cbind(t(pasture_term3_raw[,1,]), 
                       t(pasture_term3_raw[,2,]), 
                       t(pasture_term3_raw[,3,]),
                       t(pasture_term3_raw[,4,]))

# 6 species in the phylogeny represent 2 species in contemporary taxonomy (i.e. 
# map to 12 species)
# map these 906 phylogenetic species back onto the 912 species in the dataset
phylo_rnames <- gsub(".*\\[(.*)\\,pasture\\]", "\\1", row.names(pasture_term3))
phylo_indexing <- match(uc_full$phylo[uc_index], phylo_rnames)

# combine fixed effects and two random effects
pasture_term <- pasture_term1 + pasture_term2 + pasture_term3[phylo_indexing,]

## save coefficients ----
# store outputs in single object
out <- list(species = uc_full$species[uc_index],
            phylo = uc_full$phylo[uc_index],
            lpo = lpo,
            relev_term = relev_term,
            relev2_term = relev2_term,
            woody_veg_sc_term = woody_veg_sc_term, 
            pasture = pasture_term)

## pass elev_median, elev_breadth, and woody_veg scaling through here?
saveRDS(out, "outputs/lpo_and_coefs_EC.rds")

## spatially varying terms ----
# spatial_input_layers <- readRDS("outputs/spatial_input_layers.rds")
# generic_input_layers <- readRDS("outputs/generic_input_layers.rds")
# 
# spatial_input_layers[, `:=`(elev_median = NULL, 
#                             elev_breadth = NULL)]
# 
# # get scaling for woody vegetation
# wv_centering <- attr(cu_cov$woody_veg_sc, "scaled:center")
# wv_scaling <- attr(cu_cov$woody_veg_sc, "scaled:scale")
# wv_val_forest <- cu_cov %>% 
#     filter(not_forest == -1) %>%
#     pull(woody_veg_sc) %>%
#     unique %>% as.vector
# 
# generic_input_layers[, pct_cell_pasture := pct_unmasked * (1-pct_forest)]
# generic_input_layers[, pct_cell_forest := pct_unmasked * pct_forest]
# generic_input_layers[, woody_veg_sc := (avg_treecover - wv_centering)/wv_scaling]
# 
# coefs <- data.table(species = out$species, 
#                     relev_term = out$relev_term[,1], 
#                     relev2_term = out$relev2_term[,1], 
#                     woody_veg_sc_term = out$woody_veg_sc_term[,1], 
#                     lpo_forest = out$lpo_forest[,1], 
#                     lpo_pasture = out$lpo_pasture[,1])
# 
# spatial_input_layers[generic_input_layers, `:=`(
#     woody_veg_sc = i.woody_veg_sc, 
#     pct_cell_forest = i.pct_cell_forest, 
#     pct_cell_pasture = i.pct_cell_pasture
# ), on = "index"]
# 
# spatial_input_layers[coefs, `:=`(
#     logit_occ_forest = lpo_forest + woody_veg_sc_term * wv_val_forest + 
#         relev_term * relev + relev2_term * relev^2, 
#     lpo_pasture = i.lpo_forest, 
#     logit_occ_pasture = lpo_pasture + woody_veg_sc_term * woody_veg_sc  + 
#         relev_term * relev + relev2_term * relev^2), 
#     on = "species"]
# 
# any(is.na(sp_i$relev))
# # check <- readRDS("outputs/stars_lookup.rds") %>%
# #     mutate(forest_occ = NA, pasture_occ = NA)
# # sp_i <- spatial_input_layers[species == "Zonotrichia_capensis",]
# sp_i <- spatial_input_layers[species == "Adelomyia_melanogenys",]
# 
# sp_i <- spatial_input_layers %>% filter(species == "Adelomyia_melanogenys")
# sp_i <- spatial_input_layers %>% filter(species == "Aburria_aburri")
# 
# check[["forest_occ"]][sp_i$index] <- sp_i$woody_veg_sc#boot::inv.logit(sp_i$logit_occ_forest)
# check[["pasture_occ"]][sp_i$index] <- sp_i$lpo_pasture#boot::inv.logit(sp_i$logit_occ_pasture)
# # 
# # 
# plot(check["forest_occ"], breaks = "equal")
# plot(check["pasture_occ"], breaks = "equal")
# 
# spatial_input_layers[, `:=`(
#     units_forest = (pct_cell_forest * 1000^2)/(100^2*pi), 
#     units_pasture = (pct_cell_pasture * 1000^2)/(100^2*pi)
#     )]
# 
# spatial_input_layers[, `:=`(
#     amt_occ = units_forest * boot::inv.logit(logit_occ_forest) + 
#         units_pasture * boot::inv.logit(logit_occ_pasture)
# )]
# 
# spatial_input_layers[,prop_range := amt_occ/sum(amt_occ), by="species"]
# 
# new <- spatial_input_layers[, sum(prop_range), by="index"]
# new <- species_richness <- spatial_input_layers[, sum(amt_occ >= 1), by="index"]
# 
# 
# check <- readRDS("outputs/stars_lookup.rds") %>%
#     mutate(sum_range = NA)
# 
# check[["sum_range"]][new$index] <- new$V1
# check <- check %>%
#     mutate(sum_range_resc = sum_range/mean(sum_range, na.rm=T))
# check
# plot(check["sum_range"])
# plot(check["sum_range_resc"], na.rm=T)
# 
# hist(log(check[["sum_range"]]/mean(check[["sum_range"]], na.rm=T)))
# exp(2)
# 
# 
# sp_i <- spatial_input_layers %>% 
#     filter(species == "Aburria_aburri")
# 
