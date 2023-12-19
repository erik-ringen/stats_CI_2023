library(rethinking)
library(tidyverse)
data("reedfrogs")

d <- reedfrogs

# make the tank cluster variable
d$tank <- 1:nrow(d)

dat <- list(
  S = d$surv,
  N = d$density,
  tank = d$tank
)

# start with basic varying intercepts model
m1 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank],
    # priors
    a[tank] ~ dnorm(a_bar, sigma_tank),
    a_bar ~ dnorm(0, 1.5),
    sigma_tank ~ dexp(1)
  ),
  data = dat,
  chains = 4,
  cores = 4,
  log_lik = TRUE
)

# Now, try adding size and predation predictors
dat$size <- ifelse(d$size == "small", 0, 1)
dat$pred <- ifelse(d$pred == "no", 0, 1)

# size only
m2 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank] + b_s*size,
    # priors
    a[tank] ~ dnorm(a_bar, sigma_tank),
    a_bar ~ dnorm(0, 1.5),
    sigma_tank ~ dexp(1),
    b_s ~ dnorm(0, 0.5)
  ),
  data = dat,
  chains = 4,
  cores = 4,
  log_lik = TRUE
)

# predation only
m3 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank] + b_p*pred,
    # priors
    a[tank] ~ dnorm(a_bar, sigma_tank),
    a_bar ~ dnorm(0, 1.5),
    sigma_tank ~ dexp(1),
    b_p ~ dnorm(0, 0.5)
  ),
  data = dat,
  chains = 4,
  cores = 4,
  log_lik = TRUE
)

# both predictors
m4 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank] + b_s*size + b_p*pred,
    # priors
    a[tank] ~ dnorm(a_bar, sigma_tank),
    a_bar ~ dnorm(0, 1.5),
    sigma_tank ~ dexp(1),
    b_s ~ dnorm(0, 0.5),
    b_p ~ dnorm(0, 0.5)
  ),
  data = dat,
  chains = 4,
  cores = 4,
  log_lik = TRUE
)

post <- extract.samples(m4)

# interaction between predictors
m5 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank] + b_s*size + b_p*pred + b_s_p*size*pred,
    # priors
    a[tank] ~ dnorm(a_bar, sigma_tank),
    a_bar ~ dnorm(0, 1.5),
    sigma_tank ~ dexp(1),
    b_s ~ dnorm(0, 0.5),
    b_p ~ dnorm(0, 0.5),
    b_s_p ~ dnorm(0, 0.5)
  ),
  data = dat,
  chains = 4,
  cores = 4,
  log_lik = TRUE
)

# Check out how sigma_tank changes when we add more predictors
post_m1 <- extract.samples(m1)
post_m5 <- extract.samples(m5)

dens(post_m1$sigma_tank, xlim = c(0, 3), ylim = c(0, 3), col = "orange", lwd = 2, xlab = "sigma_tank")
dens(post_m5$sigma_tank, add=T, col = "skyblue", lwd = 2)

# LOOCV model comparison
compare(m1, m2, m3, m4, m5, func = PSIS)

loo_m5 <- PSIS(m5, pointwise = T)

preds <- sim(m5)
pred_mu <- apply(preds, 2, mean)
pred_PI <- apply(preds, 2, PI)

# Organize model sims into dataframe
d_sims <- data.frame(
  pred_mu = pred_mu,
  lower = pred_PI[1,],
  upper = pred_PI[2,],
  tank = dat$tank,
  density = dat$N,
  size = d$size,
  pred = d$pred,
  pareto_k = loo_m5$k,
  prop_surv = d$surv/d$density
)

ggplot(d_sims, aes(x = tank, y = pred_mu/density)) +
  facet_wrap(size ~ pred, scales = "free_x") +
  geom_point(color = ifelse(d_sims$pareto_k > 1, "red", "cornflowerblue")) +
  geom_errorbar(aes(ymin = lower/density, ymax = upper/density), width = 0, color = ifelse(d_sims$pareto_k > 1, "red", "cornflowerblue")) +
  geom_point(aes(y = prop_surv, x = tank)) +
  ylab("Proportion Survive")  

## Is there more between-tank variance when predation occurs (heteroskedasticity)?
dat$pred_id <- case_when(
  dat$pred == 0 & dat$size == 0 ~ 1,
  dat$pred == 1 & dat$size == 0 ~ 2,
  dat$pred == 0 & dat$size == 1 ~ 3,
  dat$pred == 1 & dat$size == 1 ~ 4
)

m6 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a_bar + z[tank]*sigma_tank_pred[pred_id] + b_s*size + b_p*pred + b_s_p*size*pred,
    # priors, using non-centered parameterization
    z[tank] ~ dnorm(0, 1),
    a_bar ~ dnorm(0, 1.5),
    sigma_tank_pred[pred_id] ~ dexp(1),
    b_s ~ dnorm(0, 0.5),
    b_p ~ dnorm(0, 0.5),
    b_s_p ~ dnorm(0, 0.5)
  ),
  data = dat,
  chains = 4,
  cores = 4,
  log_lik = TRUE
)

loo_m6 <- PSIS(m6, pointwise = T, n = 1e5)
sum(loo_m6$k > 1)

summary(m6)

compare(m1, m2, m3, m4, m5, m6, func = PSIS)

# Now lets calculate the causal effects of predation and tadpole size
post_m6 <- extract.samples(m6)

# function to get the probability of survival under different conditions, marginalizing over tanks
epred_m6 <- function(post, pred, size){
  prob <- with(post, inv_logit(a_bar + b_s*size + b_p*pred + b_s_p*size*pred))
               
  return(prob)   
}

prob_small_nopred <- epred_m6(post_m6, pred = 0, size = 0)
prob_small_pred <- epred_m6(post_m6, pred = 1, size = 0)
prob_big_nopred <- epred_m6(post_m6, pred = 0, size = 1)
prob_big_pred <- epred_m6(post_m6, pred = 1, size = 1)

# Causal effect of predation: depends on size of tadpoles
diff_pred_small <- prob_small_pred - prob_small_nopred
diff_pred_big <- prob_big_pred - prob_big_nopred

dens(diff_pred_small, xlim = c(-1, 1), col = "orange", xlab = "Diff (Prob Survive)")
dens(diff_pred_big, add = T, col = "indianred")

# Causal effect of tadpole size: depends on predation level
diff_size_nopred <- prob_big_nopred - prob_small_nopred
diff_size_pred <- prob_big_pred - prob_small_pred

dens(diff_size_nopred, xlim = c(-1, 1), col = "skyblue", xlab = "Diff (Prob Survive)")
dens(diff_size_pred, add = T, col = "royalblue")
