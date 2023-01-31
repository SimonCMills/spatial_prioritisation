# A minimal example to demonstrate FLM approach when effects covary along a 
# gradient 

# packages
library(dplyr); library(ggplot2); library(mgcv)

# simulate data
n_obs <- 30
n_vars <- 10

mat <- matrix(rnorm(n_obs*n_vars), ncol=n_vars)
betas <- sin(seq(-2, 2, len=n_vars)) # correlation in effects across vars
response <- mat %*% betas + rnorm(n_obs, 0, 1)

# Multiple regression ----
fit_lm <- lm(response ~ mat)

preds_lm <- summary(fit_lm)$coefficients %>%
    as.data.frame %>%
    slice(-1) %>%
    rename(fit = Estimate, 
           se = `Std. Error`) %>%
    mutate(id_var = seq_len(n_vars), 
           lwr = fit - 2*se, 
           upr = fit + 2*se)

# Multiple models ----
catch <- c()
for(i in seq_len(n_vars)) {
    catch[[i]] <- summary(lm(response ~ mat[,i]))$coefficients[2,]
}

pred_multiple <- bind_rows(catch) %>%
    rename(fit = Estimate, 
           se = `Std. Error`) %>%
    mutate(id_var = seq_len(n_vars), 
           lwr = fit - 2*se, 
           upr = fit + 2*se)

# GAM ---- 
# format data 
df <- data.frame(response)
lags <- matrix(rep(seq_len(n_vars), each=n_obs), ncol=n_vars)
df$lags <- lags
df$covar <- mat

fit_gam <- gam(response ~ s(lags, by=covar, k=5, bs="cs"), data=df)

preds_newdf <- tibble(lags = seq_len(n_vars), 
                covar = 1)

preds <- predict(fit_gam, newdata=preds_newdf, se = T, type="terms")# %>%
    
pred_gam <- data.frame(fit = preds$fit[,1], se = preds$se.fit[,1]) %>%
    mutate(id_var = seq_len(n_vars), 
           lwr = fit - 2*se, 
           upr = fit + 2*se)

# Plot ----
ggplot(pred_gam, aes(id_var-.1, fit, ymax=upr, ymin=lwr)) +
    # estimates from FLM (black)
    geom_point() + 
    geom_linerange() +
    # estimates from multiple models (blue)
    geom_point(data=pred_multiple, aes(id_var), col="blue") +
    geom_linerange(data=pred_multiple, aes(id_var), col="blue") +
    # estimates from single model with all covariates (red)
    geom_point(data=preds_lm, aes(x = id_var + .1), col="red") +
    geom_linerange(data=preds_lm, aes(x = id_var + .1), col="red") +
    geom_linerange(aes(xmin = id_var - .5, xmax = id_var + .5, y=betas), 
                   inherit.aes = FALSE, linetype = 3) +
    scale_x_continuous(breaks = seq_len(n_vars)) +
    labs(x = "Variable id",
         y = "Estimate (95% CI)",
         caption = paste(strwrap("Estimates from FLM (blue), multiple models (blue), 
         and multiple regression (red), plus truth (dashed lines)"), collapse=" "))
