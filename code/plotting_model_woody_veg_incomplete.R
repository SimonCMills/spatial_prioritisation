# initial exploration of EC model (first draft; there are issues with it)

## housekeeping ----
library(brms); library(flocker); library(dplyr)

## rebuild fit from incomplete fit ----
A <- readRDS("stan_out/EC_model_run5_woody_veg_effs/A.rds")
fd <- readRDS("stan_out/EC_model_run5_woody_veg_effs/fd.rds")

rebuilt_fit_fname <- "stan_out/EC_model_run6_woody_veg_effs/rebuilt_fit.rds"
if(!file.exists(rebuilt_fit_fname)) {
    fit_fnames <- list.files("stan_out/EC_model_run6_woody_veg_effs/", 
                             ".csv",
                             full.names = TRUE)
    
    fit_samples <- lapply(fit_fnames, brms:::read_csv_as_stanfit)
    
    # identify the minimum number of completed draws across chains
    min_draws <- min(sapply(fit_samples, function(x) x@sim$n_save))
    
    # subset each chain down to minimum number of draws
    for(i in seq_along(fit_samples)) {
        fit_samples[[i]]@sim$samples[[1]] <- fit_samples[[i]]@sim$samples[[1]][1:min_draws,]
        attributes(fit_samples[[i]]@sim$samples[[1]])$sampler_params <-
            attributes(fit_samples[[i]]@sim$samples[[1]])$sampler_params[1:min_draws,]
        fit_samples[[i]]@sim$n_save <- min_draws
        fit_samples[[i]]@sim$iter <- min_draws 
        # don't need to change this entry and will never need to refer to it, but
        # note that it will contain incorrect info after rebuild
        # fit_samples[[i]]@stan_args[[1]]$iter <- min_draws
    }
    
    # convert list of stanfits to single stanfit
    sfits <- rstan::sflist2stanfit(fit_samples)
    
    fit_empty <- flock(
        # occupancy
        f_occ = ~ 1 + forest_dep + pasture +
            pasture:forest_dep +
            woody_veg_sc + 
            woody_veg_sc:forest_dep +
            # dietary traits
            diet_inv + diet_inv:pasture + diet_inv:woody_veg_sc +
            diet_carn + diet_carn:pasture + diet_carn:woody_veg_sc +
            diet_fruitnect + diet_fruitnect:pasture + diet_fruitnect:woody_veg_sc +
            diet_gran + diet_gran:pasture + diet_gran:woody_veg_sc +
            # mass and foraging traits
            log_mass_sc + log_mass_sc:pasture + log_mass_sc:woody_veg_sc +
            strata_index + strata_index:pasture + strata_index:woody_veg_sc +
            strata_aerial + strata_aerial:pasture + strata_aerial:woody_veg_sc +
            # elevational breadth
            elev_breadth_sc + elev_breadth_sc:pasture + elev_breadth_sc:woody_veg_sc +
            # elevational position
            elev_median_sc + 
            elev_median_sc:forest_dep +
            elev_median_sc:pasture + 
            elev_median_sc:woody_veg_sc + 
            elev_median_sc:pasture:forest_dep +
            elev_median_sc:woody_veg_sc:forest_dep + 
            # elevational range
            relev + relev2 + 
            # ranefs
            (1 + pasture + woody_veg_sc + relev + relev2|species) + 
            (1|sr) + (1|sp_sr) + (1|sp_cl) +
            (1 + pasture|gr(phylo, cov = A)), # try to fit phylo woody veg too?
        # detection            
        f_det = ~ 1 + time + observer_pt + strata_index + strata_aerial + pasture + 
            woody_veg_sc + (1 + pasture|species) + (1|obs_sp), 
        # flocker info            
        data2 = list(A = A),
        flocker_data = fd, 
        rep_constant = F, 
        # run & save info 
        empty = TRUE)
    
    
    # stick the rebuilt stanfit back into brm's innards.
    fit_empty$fit <- sfits
    fit <- rename_pars(fit_empty)
    saveRDS(fit, rebuilt_fit_fname)
} else {
    fit <- readRDS(rebuilt_fit_fname)
}
# fit <- readRDS(rebuilt_fit_fname)
print(fit)

## species list for Emma ----
# lookup <- read.csv("../Colombia/colombia_birds/outputs/initial_species_list.csv")
# species_list <- data.frame(HBW = unique(pred_data$species)) %>%
#     mutate(HBW = gsub("_", " ", HBW)) %>%
#     left_join(., lookup) %>%
#     as_tibble
# saveRDS(species_list, "species_list_for_Emma.rds")

## predictions ----
pred_data_species <- fd$data %>%
  select(species, strata_index, strata_aerial, forest_dep,
         diet_inv, diet_carn, diet_fruitnect, diet_gran,
         log_mass_sc, elev_breadth, elev_breadth_sc, elev_median, elev_median_sc,
         phylo) %>%
  unique

# note: with this method only exact to 3rd decimal place
wveg_scale <- fd$data %>%
    select(point, species, woody_veg) %>%
    unique %>%
    pull(woody_veg) %>%
    sd
wveg_mean <- fd$data %>%
    select(point, species, woody_veg) %>%
    unique %>%
    pull(woody_veg) %>%
    mean

woody_veg_vals <- fd$data %>%
    select(point, woody_veg, woody_veg_sc) %>%
    unique %>%
    mutate(woody_veg_sc2 = (woody_veg - wveg_mean)/wveg_scale)

ggplot(woody_veg_vals, aes(woody_veg)) + 
    geom_histogram(boundary = 0, binwidth=.05)

ggplot(woody_veg_vals, aes(woody_veg_sc, woody_veg_sc2))
wf_vals <- c("low_wf" = (.05 - wveg_mean)/wveg_scale,
  "high_wf" = (.4 - wveg_mean)/wveg_scale,
  "forest" = (1 - wveg_mean)/wveg_scale)


woody_veg_range <- fit$data %>%
  group_by(pasture) %>%
  summarise(woody_veg_sc = quantile(woody_veg_sc, c(.1, .9))) %>%
  ungroup %>%
  slice(-1) # drop the first row so there is only one value for forest

pred_data_point <- expand.grid(species = pred_data_species$species,
                               observer_pt = "SCM",
                               elev_ALOS = seq(1000, 3500, 200),
                               pasture = c(-1, 1),
                               obs_sp = fit$data$obs_sp[1],
                               paramo = -1,
                               time = 0,
                               # relev = 0,
                               # relev2 = 0,
                               sr = fit$data$sr[1],
                               sp_sr = fit$data$sp_sr[1],
                               sp_cl = fit$data$sp_cl[1])


pred_data <- pred_data_point %>%
  full_join(., woody_veg_range) %>%
  left_join(., pred_data_species) %>%
  mutate(relev = (elev_ALOS - elev_median)/elev_breadth * 1.61,
         relev2 = relev^2)

preds <- fitted_flocker(fit, type = "occupancy", CI = c(.1, .9),
                        new_data = pred_data,
                        re_formula = ~ (1 + pasture + woody_veg_sc + relev + relev2|species) +
                          (1 + pasture|gr(phylo, cov = A)), 
                        summarise = T, ndraws = 500)
pred_full <- bind_cols(pred_data, preds) %>%
    as_tibble() #%>%
    # reshape2::melt(., id.vars = c(""))


# plot fixefs 
fixefs <- fixef(fit) %>%
    as_tibble(., rownames = "parameter") %>%
    mutate(component = c("occ", "det")[grepl("occ", parameter)+1])

ggplot(fixefs, aes(Estimate, parameter, xmin = Q2.5, xmax=Q97.5)) +
    geom_point() +
    geom_linerange() +
    geom_vline(xintercept = 0) +
    facet_wrap(~component)


# saveRDS(pred_full, "dataset_for_Zach_16-02-23.rds")

# preds ----

library(ggplot2)

preds <- fitted_flocker(fit, type = "occupancy", CI = c(.1, .9),
                        new_data = pred_data,
                        re_formula = ~ (1 + pasture + woody_veg_sc + relev + relev2|species) +
                            (1 + pasture|gr(phylo, cov = A)), 
                        summarise = FALSE, ndraws = 500)

pred_full <- bind_cols(pred_data, preds) %>%
    as_tibble() %>%
    mutate(
        wveg_round = round(woody_veg_sc, 2), 
        point_type = case_when(pasture == -1 ~ "forest", 
                               pasture == 1 & wveg_round == -1.40 ~ "low-wf pasture",
                               pasture == 1 & wveg_round == -0.67 ~ "high-wf pasture")
    ) %>%
    reshape2::melt(.,
                   measure.vars = paste0("iter_", 1:500),
                   id.vars = c("species", "elev_ALOS", "point_type", "forest_dep")) %>%
    mutate(id_iter = as.integer(as.factor(variable)), 
           id_iter_species = interaction(id_iter, species, elev_ALOS))

# n_species * n_iter
n_sp_iter <- length(levels(pred_full$id_iter_species))

pred_full2 <- pred_full %>%
    mutate(value2 = boot::inv.logit(
        boot::logit(value) + rnorm(n_sp_iter, 0, 2.7)[id_iter_species]
    )
    )

saveRDS(pred_full2, "outputs/pred_full2.rds")

library(ggplot2)
summ <- pred_full2 %>%
    group_by(elev_ALOS, point_type, variable) %>%
    summarise(SR = sum(value2)) %>%
    group_by(elev_ALOS, point_type) %>%
    summarise(estimate = mean(SR), 
              lwr = quantile(SR, .05), 
              upr = quantile(SR, .95))


ggplot(summ, aes(elev_ALOS, estimate, ymin = lwr, ymax = upr)) +
    geom_line(aes(col = point_type)) +
    geom_ribbon(aes(fill = point_type), alpha=.2)


summ2 <- pred_full2 %>%
    group_by(elev_ALOS, point_type, variable, forest_dep) %>%
    summarise(SR = sum(value2)) %>%
    group_by(elev_ALOS, point_type, forest_dep) %>%
    summarise(estimate = mean(SR), 
              lwr = quantile(SR, .05), 
              upr = quantile(SR, .95))

ggplot(summ2, aes(elev_ALOS, estimate, ymin = lwr, ymax = upr)) +
    geom_line(aes(col = point_type)) +
    geom_ribbon(aes(fill = point_type), alpha=.2) +
    facet_wrap(~forest_dep)

summ_sp <- pred_full2 %>%
    group_by(point_type, species) %>%
    summarise(n_pt = sum(value2)) %>%
    group_by(species, point_type) %>%
    summarise(estimate = mean(n_pt), 
              lwr = quantile(n_pt, .05), 
              upr = quantile(n_pt, .95))


ggplot(summ_sp, aes(point_type, estimate)) + geom_violin() +
    scale_y_continuous(trans="log", breaks = 2^(0:10))


## sparing-sharing ----

woody_veg_range2 <- woody_veg_range %>%
    mutate(woody_veg_sc = wf_vals[3:1]) %>%
    bind_rows(., tibble(pasture = 1, woody_veg_sc = (0-wveg_mean)/wveg_scale))

pred_data_wf <- pred_data_point %>%
    full_join(., woody_veg_range2) %>%
    left_join(., pred_data_species) %>%
    mutate(relev = (elev_ALOS - elev_median)/elev_breadth * 1.61,
           relev2 = relev^2)

preds_wf <- fitted_flocker(fit, type = "occupancy", CI = c(.1, .9),
                        new_data = pred_data_wf,
                        re_formula = ~ (1 + pasture + woody_veg_sc + relev + relev2|species) +
                            (1 + pasture|gr(phylo, cov = A)), 
                        summarise = FALSE, ndraws = 500)

pred_full_wf <- bind_cols(pred_data_wf, preds_wf) %>%
    as_tibble() %>%
    mutate(
        wveg_round = round(woody_veg_sc, 2), 
        point_type = case_when(pasture == -1 ~ "forest", 
                               pasture == 1 & wveg_round < -1 & wveg_round > -1.4 ~ "low_wf_pasture",
                               pasture == 1 & wveg_round < -1.4 ~ "zero_wf_pasture",
                               pasture == 1 & wveg_round > -1 ~ "high_wf_pasture")
    ) %>%
    reshape2::melt(.,
                   measure.vars = paste0("iter_", 1:500),
                   id.vars = c("species", "elev_ALOS", "point_type", "forest_dep")) %>%
    mutate(id_iter = as.integer(as.factor(variable)), 
           id_iter_species = interaction(id_iter, species, elev_ALOS))

# n_species * n_iter
n_sp_iter <- length(levels(pred_full$id_iter_species))

pred_full2_wf <- pred_full_wf %>%
    mutate(value2 = boot::inv.logit(
        boot::logit(value) + rnorm(n_sp_iter, 0, 2.7)[id_iter_species]
    )
    )

saveRDS(pred_full2_wf, "outputs/pred_full2_wf.rds")

# SR
n_pt_lscape <- 50
sparing_sharing <- pred_full2_wf[pred_full2_wf$point_type == "forest",]
sparing_sharing$sparing_high <- 
    1 - 
    (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "forest"])^(.05 * n_pt_lscape) * # probability of not being on 1 point
    (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "zero_wf_pasture"])^(.95 * n_pt_lscape) # p(not on 19 points)

sparing_sharing$sharing_high <- 1 - (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "low_wf_pasture"])^(n_pt_lscape)

sparing_sharing$sparing_low <- 1 - 
    (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "forest"])^(.4 * n_pt_lscape) * 
    (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "zero_wf_pasture"])^(.6 * n_pt_lscape)

sparing_sharing$sharing_low <- 1 - (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "high_wf_pasture"])^n_pt_lscape

# N PT
# sparing_sharing <- pred_full2_wf[pred_full2_wf$point_type == "forest",]
# sparing_sharing$sparing_high <- 
#     pred_full2_wf$value2[pred_full2_wf$point_type == "forest"] + # probability of not being on 1 point
#     pred_full2_wf$value2[pred_full2_wf$point_type == "zero_wf_pasture"] # p(not on 19 points)
# 
# sparing_sharing$sharing_high <- 1 - (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "low_wf_pasture"])^20
# 
# sparing_sharing$sparing_low <- 1 - 
#     (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "forest"])^8 * 
#     (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "zero_wf_pasture"])^12
# 
# sparing_sharing$sharing_low <- 1 - (1 - pred_full2_wf$value2[pred_full2_wf$point_type == "high_wf_pasture"])^20



summ_high <- sparing_sharing %>%
    group_by(elev_ALOS, point_type, variable) %>%
    summarise(SR_spr = sum(sparing_high), 
              SR_shr = sum(sharing_high),
              SR_diff = SR_spr - SR_shr) %>%
    group_by(elev_ALOS, point_type) %>%
    summarise_at(vars(SR_spr, SR_shr, SR_diff), 
                 .funs = list(mean = mean, 
                              lwr = function(x) quantile(x, .05), 
                              upr = function(x) quantile(x, .95))
    )

summ_low <- sparing_sharing %>%
    group_by(elev_ALOS, point_type, variable) %>%
    summarise(SR_spr = sum(sparing_low), 
              SR_shr = sum(sharing_low),
              SR_diff = SR_spr - SR_shr) %>%
    group_by(elev_ALOS, point_type) %>%
    summarise_at(vars(SR_spr, SR_shr, SR_diff), 
                 .funs = list(mean = mean, 
                              lwr = function(x) quantile(x, .05), 
                              upr = function(x) quantile(x, .95))
    )

# repackage
new_colnames <- c("elev", "mean", "lwr", "upr")
n_elevs <- length(unique(summ_high$elev_ALOS))

summ2_high <- bind_rows(select(summ_high, contains("spr")) %>% setNames(new_colnames), 
                        select(summ_high, contains("shr")) %>% setNames(new_colnames), 
                        select(summ_high, contains("diff")) %>% setNames(new_colnames)) %>%
    cbind(type = rep(c("spr", "shr", "diff"), each = n_elevs))

summ2_low <- bind_rows(select(summ_low, contains("spr")) %>% setNames(new_colnames), 
                       select(summ_low, contains("shr")) %>% setNames(new_colnames), 
                       select(summ_low, contains("diff")) %>% setNames(new_colnames)) %>%
    cbind(type = rep(c("spr", "shr", "diff"), each = n_elevs))


p1 <- summ2_high %>%
    filter(type != "diff") %>%
    ggplot(aes(elev, mean, ymin=lwr, ymax=upr)) +
    geom_line(aes(col=type)) +
    geom_ribbon(aes(fill = type), alpha=.2)

p2 <- summ2_high %>%
    filter(type == "diff") %>%
    ggplot(aes(elev, mean, ymin=lwr, ymax=upr)) +
    geom_line(aes(col=type)) +
    geom_ribbon(aes(fill = type), alpha=.2) +
    geom_hline(yintercept = 0, lty = "longdash")

p_both <- egg::ggarrange(p1, p2)
# dir.create("figures")
ggsave("figures/sparing_sharing_50_pts_high.png", plot = p_both)

p1 <- summ2_low %>%
    filter(type != "diff") %>%
    ggplot(aes(elev, mean, ymin=lwr, ymax=upr)) +
    geom_line(aes(col=type)) +
    geom_ribbon(aes(fill = type), alpha=.2)

p2 <- summ2_low %>%
    filter(type == "diff") %>%
    ggplot(aes(elev, mean, ymin=lwr, ymax=upr)) +
    geom_line(aes(col=type)) +
    geom_ribbon(aes(fill = type), alpha=.2) +
    geom_hline(yintercept = 0, lty = "longdash")

egg::ggarrange(p1, p2)
