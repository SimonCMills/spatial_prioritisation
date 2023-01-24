fit <- readRDS("warmup_27-10-2022_1032_j2789255.rds")

fit2 <- fit
thinning <- 5
for(i in seq_len(3)) {
  fit2$fit@sim$samples[[i]] <- fit$fit@sim$samples[[i]][seq(1, 500, thinning),]
  attr(fit2$fit@sim$samples[[i]], "sampler_params") <-
    attr(fit$fit@sim$samples[[i]], "sampler_params")[seq(1, 500, thinning),]
}

# fit2$fit@sim$warmup2 <- rep(0, 4)
fit2$fit@sim$thin <- thinning

fit2$fit@sim$n_save <- rep(100, 3)

ndraws(fit)
ndraws(fit2)

attr(fit$fit@sim$samples[[1]], "sampler_params")
attr(fit2$fit@sim$samples[[1]], "sampler_params")

attr(fit2$fit@sim$samples[[1]], "sampler_params")
fit2$fit@sim$samples[[1]]

df <- fit2$data

names(df)

pred_data_species <- cu_cov %>%
  select(species, strata_index, strata_aerial, forest_dep,
         diet_inv, diet_carn, diet_fruitnect, diet_gran,
         log_mass_sc, elev_breadth, elev_breadth_sc, elev_median, elev_median_sc,
         phylo) %>%
  unique

woody_veg_range <- df %>%
  # filter(not_forest == 1) %>%
  group_by(not_forest) %>%
  summarise(woody_veg_sc = quantile(woody_veg_sc, c(.1, .9))) %>%
  ungroup %>%
  slice(-1)

pred_data_point <- expand.grid(species = pred_data_species$species,
                               observer_pt = "SCM",
                               elev_ALOS = seq(500, 3500, 500),
                               not_forest = c(-1, 1),
                               obs_sp = df$obs_sp[1],
                               paramo = -1,
                               time = 0,
                               # relev = 0,
                               # relev2 = 0,
                               sr = df$sr[1],
                               sp_sr = df$sp_sr[1],
                               sp_cl = df$sp_cl[1])


pred_data <- pred_data_point %>%
  full_join(., woody_veg_range) %>%
  left_join(., pred_data_species) %>%
  mutate(relev = (elev_ALOS - elev_median)/elev_breadth * 1.61,
         relev2 = relev^2)

library(flocker)

preds <- fitted_flocker(fit2, type = "occupancy", new_data = pred_data,
                        re_formula = ~ (1 + not_forest + paramo + woody_veg_sc + relev + relev2|species) +
                          (1 + not_forest|gr(phylo, cov = A)))
strata_index

pred_full <- bind_cols(pred_data, preds)

cu_cov <- readRDS("cu_cov.rds")
temp <- readRDS("temp_file.rds")

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
  group_by(elev_ALOS, not_forest, woody_veg_sc) %>%
  summarise(SR = sum(estimate))


summ2 <- pred_full %>%
  # filter(elev_median < 3000) %>%
  group_by(species, forest_dep, elev_median, elev_ALOS) %>%
  summarise(RR1 = (estimate[not_forest == -1]/estimate[not_forest == 1][1]),
            RR2 = (estimate[not_forest == -1]/estimate[not_forest == 1][2])) %>%
  group_by(species, forest_dep, elev_median) %>%
  summarise(RR1 = mean(RR1),
            RR2 = mean(RR2))

ggplot(summ, aes(elev_ALOS, SR, col=not_forest)) + geom_point()

## RR1 is the comparison with wildlife unfriendly pasture
summ2 %>%
  left_join(., temp) %>%
ggplot(aes(elev_median, (RR1))) + geom_point(aes(col=in_range < 0 )) + facet_wrap(~forest_dep) +
  scale_y_continuous(trans = "log", breaks = 10^(-5:10)) + geom_hline(yintercept = 1) +
  stat_smooth() +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

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
