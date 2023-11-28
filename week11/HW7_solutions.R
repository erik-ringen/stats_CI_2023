library(rethinking)
library(dagitty)

## 8H4
data(nettle)
d <- nettle
d$lang.per.cap <- d$num.lang/d$k.pop
d$log_lang_per_cap <- log(d$lang.per.cap)
d$log_area <- log(d$area)

### Drawing the Bayesian Owl as an example for final project: you did not need to go into this much detail for your homework!!

## Step (1): Define theoretical estimands

# (i) The mean difference in log(lang per capita) given a +1 sd increase in the avg length of growing season.

# (ii) The mean difference in log(land per capita) given a +1 sd increase in the std dev. of growing season 

## Step (2): Scientific causal models
g <- dagitty("dag{
             mu_grow_season -> log_lang_cap; 
             log_area -> sd_grow_season -> log_lang_cap;
             log_area -> log_lang_cap}")

plot(g)
impliedConditionalIndependencies(g)

## Step (3): Build statistical models, starting with prior predictive simulations
# Start by standardizing continuous variables (z-scoring)
d$log_lang_per_cap_s <- standardize(d$log_lang_per_cap)
d$mean_growing_season_s <- standardize(d$mean.growing.season)
d$sd_growing_season_s <- standardize(d$sd.growing.season)
d$log_area_s <- standardize(d$log_area)

sim_prior_main <- function(n_obs){
  n = 1
  ### Priors ###
  # Intercepts 
  a_area = rnorm(n, mean = 0, sd = 1)
  a_mgs = rnorm(n, mean = 0, sd = 1)
  a_sdg = rnorm(n, mean = 0, sd = 1)
  a_lpc = rnorm(n, mean = 0, sd = 1)
  
  # Slopes
  b_area_lpc = rnorm(n, mean = 0, sd = 1)
  b_area_sdg = rnorm(n, mean = 0, sd = 1)
  b_mgs_lpc = rnorm(n, mean = 0, sd = 1)
  b_sdg_lpc = rnorm(n, mean = 0, sd = 1)
  
  # Residual standard deviations
  sigma_area = rexp(n, 1)
  sigma_lpc = rexp(n, 1)
  sigma_sdg = rexp(n, 1)
  sigma_mgs = rexp(n, 1)
  
  # Simulate variables
  log_area_s = rnorm(n_obs, a_area, sigma_area)
  
  mean_growing_season_s = rnorm(n_obs, a_mgs, sigma_mgs)
  
  # note how the simulated values from log_area_s are part of the linear model here
  sd_growing_season_s = rnorm(n_obs, a_sdg + b_area_sdg*log_area_s, sigma_sdg)
  
  log_lang_per_cap_s = rnorm(n_obs, a_lpc + b_area_lpc*log_area_s + b_mgs_lpc*mean_growing_season_s + b_sdg_lpc*sd_growing_season_s)
  
  # organize simulated data into a dataframe
  d_sim <- data.frame(
    log_area_s = log_area_s,
    mean_growing_season_s = mean_growing_season_s,
    sd_growing_season_s = sd_growing_season_s,
    log_lang_per_cap_s = log_lang_per_cap_s
  )
  
  # organize parameter values (which in this case correspond to our estimands), so that we can check whether we recover them during model fitting
  d_pars <- list(b_mgs_lpc = b_mgs_lpc, b_sdg_lpc = b_sdg_lpc)
  
  return(
    list(d_sim = d_sim, d_pars = d_pars)
  )
}
  
# Step 4: Simulate from statistical model to validate that we can recover our estimand, in principle
set.seed(24)
sim_dataset <- sim_prior_main(n_obs = nrow(d))

head(sim_dataset$d_sim)
sim_dataset$d_pars

# Check out the distribution of simulated data vs observed
dens(sim_dataset$d_sim$log_area_s, col = "skyblue", xlim = c(-10, 10))
dens(d$log_area_s, add = T)

dens(sim_dataset$d_sim$mean_growing_season_s, col = "skyblue", xlim = c(-10, 10))
dens(d$mean_growing_season_s, add = T)

dens(sim_dataset$d_sim$log_lang_per_cap_s, col = "skyblue", xlim = c(-10, 10))
dens(d$log_lang_per_cap_s, add = T)

# At this stage, you may refine your prior model if the predictions are really far off from the observed data.

# Try fitting the simulated data to the same model
m_prior_main <- quap(
  alist(
    # model all variables at once to simulate
    log_area_s ~ dnorm(mu_area, sigma_area),
    log_lang_per_cap_s ~ dnorm(mu_lpc, sigma_lpc),
    sd_growing_season_s ~ dnorm(mu_sdg, sigma_sdg),
    mean_growing_season_s ~ dnorm(mu_mgs, sigma_mgs),
    
    # Linear models
    mu_area <- a_area,
    mu_mgs <- a_mgs,
    mu_sdg <- a_sdg + b_area_sdg*log_area_s,
    mu_lpc <- a_lpc + b_area_lpc*log_area_s + b_mgs_lpc*mean_growing_season_s + b_sdg_lpc*sd_growing_season_s,
  
  # Priors
  c(a_area, a_mgs, a_sdg, a_lpc) ~ dnorm(0, 1),
  c(b_area_sdg, b_area_lpc, b_mgs_lpc, b_sdg_lpc) ~ dnorm(0, 1), 
  c(sigma_area, sigma_mgs, sigma_sdg, sigma_lpc) ~ dexp(1)
  ),
  # note that we are passing simulated data here
  data = sim_dataset$d_sim
)

post_sim <- extract.samples(m_prior_main)

# Compare posterior estimates with known true parameter values
dens(post_sim$b_mgs_lpc)
abline(v = sim_dataset$d_pars$b_mgs_lpc, lty = "dashed")

dens(post_sim$b_sdg_lpc)
abline(v = sim_dataset$d_pars$b_sdg_lpc, lty = "dashed")

# In research we would do more than just 1 simulation, but the principle is important: it shows that we can use a statistical model to recover the estimand. Our analysis is justified, conditional on our assumptions.

## Step (5): Analyze real data
# starting with the model we simulated from, main effects only
m_main <- quap(
  alist(
    log_lang_per_cap_s ~ dnorm(mu_lpc, sigma_lpc),
    
    # Linear model
    mu_lpc <- a_lpc + b_area_lpc*log_area_s + b_mgs_lpc*mean_growing_season_s + b_sdg_lpc*sd_growing_season_s,
    
    # Priors
    a_lpc  ~ dnorm(0, 1),
    c(b_area_lpc, b_mgs_lpc, b_sdg_lpc) ~ dnorm(0, 1), 
    sigma_lpc ~ dexp(1)
  ),
  # now we use the real data
  data = d
)

# Evaluate results
post_main <- extract.samples(m_main)
dens(post_main$b_mgs_lpc)
dens(post_main$b_sdg_lpc)

# Consider a model with an interaction. You could also do prior predictive simulations for this.
m_int <- quap(
  alist(
    log_lang_per_cap_s ~ dnorm(mu_lpc, sigma_lpc),
    
    # Linear model
    mu_lpc <- a_lpc + b_area_lpc*log_area_s + b_mgs_lpc*mean_growing_season_s + b_sdg_lpc*sd_growing_season_s + b_int*mean_growing_season_s*sd_growing_season_s,
    
    # Priors
    a_lpc  ~ dnorm(0, 1),
    c(b_area_lpc, b_mgs_lpc, b_sdg_lpc, b_int) ~ dnorm(0, 1), 
    sigma_lpc ~ dexp(1)
  ),
  # now we use the real data
  data = d
)

# Evaluate results
post_int <- extract.samples(m_interaction)
dens(post_int$b_mgs_lpc)
dens(post_int$b_sdg_lpc)
dens(post_int$b_int)

# Visualize the interaction effect
newdata <- expand.grid(
  log_area_s = 0,
  mean_growing_season_s = c(-2, 0, 2),
  sd_growing_season_s = seq(from = -2, to = 2, length.out = 30)
)
head(newdata)

mu_int <- link(m_int, data = newdata)
mean_pred <- apply(mu_int, 2, mean)
PI_pred <- apply(mu_int, 2, PI)

par(mfrow = c(1, 3))

# set y-axis limits
y_lims <- range(mu_int)

# short growing season
plot(y = mean_pred[which(newdata$mean_growing_season_s == -2)], x =  newdata$sd_growing_season_s[which(newdata$mean_growing_season_s == -2)], type = "l", xlab = "sd_growing_season_s", ylab = "log_lang_per_cap_s", ylim = y_lims)

shade(PI_pred[,which(newdata$mean_growing_season_s == -2)],  newdata$sd_growing_season_s[which(newdata$mean_growing_season_s == -2)])
mtext("mean_growing_season_s = -2")

# average growing season
plot(y = mean_pred[which(newdata$mean_growing_season_s == 0)], x =  newdata$sd_growing_season_s[which(newdata$mean_growing_season_s == 0)], type = "l", xlab = "sd_growing_season_s", ylab = "log_lang_per_cap_s", ylim = y_lims)

shade(PI_pred[,which(newdata$mean_growing_season_s == 0)],  newdata$sd_growing_season_s[which(newdata$mean_growing_season_s == 0)])
mtext("mean_growing_season_s = 0")

# long growing season
plot(y = mean_pred[which(newdata$mean_growing_season_s == 2)], x =  newdata$sd_growing_season_s[which(newdata$mean_growing_season_s == 2)], type = "l", xlab = "sd_growing_season_s", ylab = "log_lang_per_cap_s", ylim = y_lims)

shade(PI_pred[,which(newdata$mean_growing_season_s == 2)],  newdata$sd_growing_season_s[which(newdata$mean_growing_season_s == 2)])
mtext("mean_growing_season_s = +2")

# Model comparison
mod_compare <- compare(m_main, m_int, func = PSIS, n = 1e5)
mod_compare

# A possible mis-specification of the interaction: it does not appear to be symmetrical!
with(d[d$mean_growing_season_s > 0,], plot(log_lang_per_cap_s ~ sd_growing_season_s))
mtext("mean_growing_seasons_s > 0")

with(d[d$mean_growing_season_s < 0,], plot(log_lang_per_cap_s ~ sd_growing_season_s))
mtext("mean_growing_seasons_s < 0")

# Try a different parameterization (mean split for the interaction)
d$mean_growing_high <- ifelse(d$mean_growing_season_s > 0, 1, 0)

m_int_2 <- quap(
  alist(
    log_lang_per_cap_s ~ dnorm(mu_lpc, sigma_lpc),
    
    # Linear model
    mu_lpc <- a_lpc + b_area_lpc*log_area_s + b_mgs_lpc*mean_growing_season_s + b_sdg_lpc*sd_growing_season_s + b_int*mean_growing_high*sd_growing_season_s,
    
    # Priors
    a_lpc  ~ dnorm(0, 1),
    c(b_area_lpc, b_mgs_lpc, b_sdg_lpc, b_int) ~ dnorm(0, 1), 
    sigma_lpc ~ dexp(1)
  ),
  # now we use the real data
  data = d
)

precis(m_int_2)

mod_compare <- compare(m_main, m_int, m_int_2, func = PSIS, n = 1e5)
mod_compare

# Now, let's use the ensemble of models to get our final estimates
post_main <- extract.samples(m_main)
post_int <- extract.samples(m_int)
post_int_2 <- extract.samples(m_int_2)

weights <- mod_compare$weight[match(rownames(mod_compare), c("m_main", "m_int", "m_int_2"))]

b_sdg <- post_main$b_sdg_lpc * weights[1] + post_int$b_sdg_lpc * weights[2] + post_int_2$b_sdg_lpc * weights[3]

b_mgs <- post_main$b_mgs_lpc * weights[1] + post_int$b_mgs_lpc * weights[2] + post_int_2$b_mgs_lpc * weights[3]

mean(b_sdg)
HPDI(b_sdg, prob = 0.9)
mean(b_sdg < 0) # posterior probability of directional hypothesis

mean(b_mgs)
HPDI(b_mgs)
mean(b_mgs > 0) # posterior probability of directional hypothesis
  
# We could also check the implied conditional independencies
plot(g)
impliedConditionalIndependencies(g)

plot(d$mean_growing_season_s, d$sd_growing_season_s)
cor(d$mean_growing_season_s, d$sd_growing_season_s)

plot(d$log_area_s, d$mean_growing_season_s)
cor(d$log_area_s, d$mean_growing_season_s) # why are these correlated? not consistent with the implications of the DAG
