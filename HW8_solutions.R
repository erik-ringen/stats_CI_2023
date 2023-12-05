library(rethinking)
data(rugged)

### 9M1 ###
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[complete.cases(d$rgdppc_2000),]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse(dd$cont_africa == 1,1,2)

dat_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer(dd$cid)
)

m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid] + b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ),
  data = dat_slim, chains = 4, cores = 4
)

# Now with uniform prior
m9.1_uni <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid] + b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dunif(0, 1)
  ),
  data = dat_slim, chains = 4, cores = 4
)

post_9.1 <- extract.samples(m9.1)
post_9.1_uni <- extract.samples(m9.1_uni)

# compare posterior distributions of sigma
dens(post_9.1$sigma)
dens(post_9.1_uni$sigma, add = T, col = "darkred")

# compare prior distributions of sigma
dens(rexp(1e4, rate = 1), ylim = c(0,1.1), lty = "dotted", xlim = c(0, 5))
dens(runif(1e4, 0, 1), add = T, col = "darkred", lty = "dotted")
abline(v = mean(post_9.1$sigma)) # around the posterior mean of sigma, both priors are flat

# Zoom in on prior range to clarify
dens(rexp(1e5, rate = 1), ylim = c(0,1.1), lty = "dotted", xlim = c(0, max(post_9.1$sigma)))
dens(runif(1e5, 0, 1), add = T, col = "darkred", lty = "dotted")

### 9M2 ###
# Now change the prior for b[cid]
# Now with uniform prior
m9.1_b <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid] + b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dexp(0.3),
    sigma ~ dunif(0, 1)
  ),
  data = dat_slim, chains = 4, cores = 4
)

precis(m9.1, depth = 2, pars = "b")
precis(m9.1_b, depth = 2, pars = "b")

### 9M3 ###
# Trying same model with different number of warmup iterations
m9.1_10 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid] + b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ),
  data = dat_slim, chains = 4, cores = 4,
  iter = 1000,
  warmup = 10
)

traceplot_ulam(m9.1_10)

m9.1_100 <- ulam(
  alist(
    log_gdp_std ~ dnorm(mu,sigma),
    mu <- a[cid] + b[cid]*(rugged_std-0.215),
    a[cid] ~ dnorm(1,0.1),
    b[cid] ~ dnorm(0,0.3),
    sigma ~ dexp(1)
  ),
  data = dat_slim, chains = 4, cores = 4,
  iter = 1000,
  warmup = 100
)

traceplot_ulam(m9.1_100) # we see evidence of convergence by 100 warmup iterations, thus the default of iter/2 for warmup may be inefficient


### BONUS: ESS for statistics other than the central tendency ###
# install.packages("posterior")
library(posterior)

# To calculate ESS, we need to seperate the samples from each chain. rethinking combines them into a single vector
post <- extract.samples(m9.1_100)

n_samps <- length(post$sigma)
n_chains <- 4

sigma_chain_split <- matrix(post$sigma, nrow = n_samps/n_chains, ncol = n_chains)

ess_bulk(sigma_chain_split) # this is what we get from rethinking precis() function, the approximate ESS for mean, median, mode
ess_tail(sigma_chain_split) # ess for estimating the 90% CI

mcse_mean(sigma_chain_split)
mcse_median(sigma_chain_split)
mcse_quantile(sigma_chain_split) # ~3 times higher than the MCSE of the central tendency
