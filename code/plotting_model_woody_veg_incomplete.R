# initial exploration of EC model (first draft; there are issues with it)

## housekeeping ----
library(cmdstanr); library(brms); library(flocker); library(dplyr)

## rebuild fit from incomplete fit ----
rebuilt_fit_fname <- "stan_out/EC_model_run5_woody_veg_effs/rebuilt_fit.rds"
if(!file.exists(rebuilt_fit_fname)) {
    fit_fnames <- list.files("stan_out/EC_model_run5_woody_veg_effs/", 
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
    
    # build object
    A <- readRDS("stan_out/EC_model_run5_woody_veg_effs/A.rds")
    fd <- readRDS("stan_out/EC_model_run5_woody_veg_effs/fd.rds")
    
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
            (1 + pasture + paramo + woody_veg_sc + relev + relev2|species) + 
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
}
fit <- readRDS(rebuilt_fit_fname)
print(fit)

## predictions ----
pred_data_species <- fd$data %>%
  select(species, strata_index, strata_aerial, forest_dep,
         diet_inv, diet_carn, diet_fruitnect, diet_gran,
         log_mass_sc, elev_breadth, elev_breadth_sc, elev_median, elev_median_sc,
         phylo) %>%
  unique

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
                        re_formula = ~ (1 + not_forest + woody_veg_sc + relev + relev2|species) +
                          (1 + pasture|gr(phylo, cov = A)), 
                        summarise = FALSE, ndraws = 500)
pred_full <- bind_cols(pred_data, preds) %>%
    as_tibble()

saveRDS(pred_full, "dataset_for_Zach_16-02-23.rds")

## junk plotting ----
# cu_cov <- readRDS("cu_cov.rds")
# temp <- readRDS("temp_file.rds")

cu_cov %>% head

tt <- cu_cov %>%
  slice(1:1000) %>%
  ungroup %>%
  mutate(relev_new = (elev_ALOS - elev_median)/elev_breadth * 1.61)

ggplot(tt, aes(relev, relev_new)) + geom_point() + geom_abline()

summary(lm(relev~ relev_new, tt))
pnorm(1.614/2)

library(ggplot2)


summ <- pred_full %>%
  group_by(elev_ALOS, pasture, woody_veg_sc) %>%
  summarise(SR = sum(estimate))


summ2 <- pred_full %>%
  # filter(elev_median < 3000) %>%
  group_by(species, forest_dep, elev_median, elev_ALOS) %>%
  summarise(RR1 = (estimate[pasture == -1]/estimate[pasture == 1][1]),
            RR2 = (estimate[pasture == -1]/estimate[pasture == 1][2]), 
            RR3 = (estimate[pasture == 1][1]/estimate[pasture == 1][2])) %>%
  group_by(species, forest_dep, elev_median) %>%
  summarise(RR1 = mean(RR1),
            RR2 = mean(RR2), 
            RR3 = mean(RR3))

ggplot(summ, aes(elev_ALOS, SR, col=factor(pasture))) + geom_point()

## RR1 is the comparison with wildlife unfriendly pasture
summ2 %>%
  #left_join(., temp) %>%
ggplot(aes(elev_median, (RR1))) + 
    geom_point() + 
    facet_wrap(~forest_dep) +
  scale_y_continuous(trans = "log", breaks = 10^(-5:10)) +
    geom_hline(yintercept = 1) +
  stat_smooth() +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

summ2 %>%
    #left_join(., temp) %>%
    ggplot(aes(elev_median, (RR3))) + 
    geom_point() + 
    facet_wrap(~forest_dep) +
    scale_y_continuous(trans = "log", breaks = 10^(-5:10)) +
    geom_hline(yintercept = 1) +
    stat_smooth() +
    theme_bw() +
    theme(panel.grid.minor = element_blank())


pred_full %>%
    filter(species == "Zonotrichia_capensis") %>%
    ggplot(aes(relev, estimate, col=pasture, ymin=Q10, ymax=Q90)) + 
    geom_point() +
    geom_linerange()

pred_full %>%
    filter(species == "Aburria_aburri") %>%
    ggplot(aes(relev, estimate, col=pasture, ymin=Q10, ymax=Q90)) + 
    geom_point() +
    geom_linerange()


nrow(df)
nrow(cu_cov)
library(dplyr)
df %>%
  filter(abs(relev) < 0.01) %>%
  select(relev, relev2) %>% plot

df %>%
  filter(abs(relev) < 0.01) %>%
  select(relev) %>%
  mutate(relev2 = relev^2) %>% points(., col="red", cex=.5)
