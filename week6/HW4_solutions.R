library(rethinking)

### 5H1 ####
data("WaffleDivorce")

# Are divorce rates independent of marriage rates, conditional on age at marriage?
d <- data.frame(
  M = standardize(WaffleDivorce$Marriage),
  D = standardize(WaffleDivorce$Divorce),
  A = standardize(WaffleDivorce$MedianAgeMarriage)
)

m1 <- quap(
  alist(
    D ~ dnorm(mu,sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 2),
    c(bM, bA) ~ dnorm(0, 1),
    sigma ~ dexp(0.5)),
  data = d
)

plot(summary(m1))

### 5H2 
m2 <- quap(
  alist(
    # likelihood for age at first marriage
    A ~ dnorm(mu_A, sigma_A),
    # likelihood for divorce rates
    D ~ dnorm(mu_D,sigma_D),
    
    # linear model for age at first marriage
    mu_A <- a_A + bM_A*M,
    # linear model for divorce rates
    mu_D <- a_D + bA_D*A,
    
    # priors
    c(a_A, a_D) ~ dnorm(0, 2),
    c(bM_A, bA_D) ~ dnorm(0, 1),
    c(sigma_A, sigma_D) ~ dexp(0.5)),
  
  data = d
)

# Counterfactual effect of doubling marriage rate
post <- extract.samples(m2)
n_samps <- nrow(post)

# Function to simulate divorce rates, conditional on our model and different marriage rates
sim_divorce <- function(post, marriage_rate) {
  
  # standardize the value
  M <- (marriage_rate - mean(WaffleDivorce$Marriage)) / sd(WaffleDivorce$Marriage)
  
  # sim age at first marriage
  A_sim <- rnorm(n_samps, post$a_A + post$bM_A*M, post$sigma_A)
  # sim divorce rates, using simulations from previous step as input
  D_sim <- rnorm(n_samps, post$a_D + post$bM_A*A_sim, post$sigma_D)
  
  # Put predictions back onto original measurement scale
  return( (D_sim * sd(WaffleDivorce$Divorce)) + mean(WaffleDivorce$Divorce)) 
}

# Let's get Washington's marriage rate
WA_obs <- WaffleDivorce$Marriage[WaffleDivorce$Location == "Washington"]
WA_cf <- WA_obs * 2 # double of observed rate

pred_obs <- sim_divorce(post, marriage_rate = WA_obs)
pred_cf <- sim_divorce(post, marriage_rate = WA_cf)

# How would the divorce rate have been different, had the marriage rate been twice as high?
par(mfrow = c(3,1))

dens(pred_obs, xlab = paste("Divorce Rate | Marriage Rate =", WA_obs))
dens(pred_cf, xlab = paste("Divorce Rate | Marriage Rate =", WA_cf))
dens(pred_cf - pred_obs, xlab = paste("Divorce Rate | Marriage Rate =", WA_cf, "-", paste("Divorce Rate | Marriage Rate =", WA_obs)))
abline(v = 0, lty = "dashed", col = "darkred")

### 5H3 ###
data(milk)

# Assume that:
# M -> N -> K, M -> K

d <- milk
d$K <- standardize(d$kcal.per.g)
d$N <- standardize(d$neocortex.perc)
d$M <- standardize(log(d$mass))

dcc <- d[complete.cases(d$K,d$N,d$M),]

m5.7 <- quap(
  alist(
    # Likelihood for N
    N ~ dnorm(mu_N, sigma_N),
    # Likelihood for K
    K ~ dnorm(mu_K,sigma_K),
    
    # Linear model for N
    mu_N <- a_N + bM_N*M,
    # Linear model for K
    mu_K <- a_K + bM_K*M + bN_K*N,
    
    # priors
    c(a_N, a_K) ~ dnorm(0, 0.2),
    c(bM_N, bM_K, bN_K) ~ dnorm(0,0.5),
    c(sigma_K, sigma_N) ~ dexp(1)
  ), data = dcc)

# Use the same template as the previous question to compute counterfactuals. The only difference is now there are two causes of the outcome (K) to include in the simulations.

# 5H4 is very open ended, I do not provide a particular solution but you should have drawn a DAG or a few DAGs and listed their testable implications, as well as fit at least one model to the data. This is also a very culture-bound question so apologies if you are not very familiar with the culture of the U.S. South it asks about.



