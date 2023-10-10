library(rethinking)

# Set up posterior samples for problems 
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(6, size = 9, prob = p_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

set.seed(100)
samples <- sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)

# 3E1. How much posterior probability lies below p = 0.2?
mean(samples < 0.2) 
sum(samples < 0.2) / length(samples) # these two are equivalent

# 3E2. How much posterior probability lies above p = 0.8?
mean(samples > 0.8)

# 3E3. How much posterior probability lies above p = 0.8?
mean(samples > 0.2 & samples < 0.8)

# 3E4. 20% of the posterior probability lies below which value of p?
quantile(samples, probs = 0.2)

# What exactly does this mean?
sample_quantiles <- rank(samples)/length(samples)
plot(sample_quantiles ~ samples)
abline(h = 0.2, lty = "dashed", col = "darkred")
abline(v = quantile(samples, probs = 0.2), lty = "dashed", col = "darkred")

# 3E5. 20% of the posterior probability lies above which value of p?
quantile(samples, probs = 1 - 0.2)

abline(h = 0.8, lty = "dashed", col = "cornflowerblue")
abline(v = quantile(samples, probs = 1 - 0.2), lty = "dashed", col = "cornflowerblue")

# 3E6. Which values of p contain the narrowest interval equal to 66% of the posterior probability?
HPDI(samples, prob = 0.66)

# 3E7. Which values of p contain 66% of the posterior probability, assuming equal posterior probability both below and above the interval?
PI(samples, prob = 0.66) # PI function is just a convenient quantile()
quantile(samples, probs = c( (1 - 0.66)/2, 1 - (1 - 0.66)/2) )

# 3M1.
likelihood <- dbinom(8, size = 15, prob = p_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

# 3M2. Draw 10,000 samples from the grid approximation from above.Then use the samples to calculate the 90% HPDI for p.
samples <- sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)
HPDI(samples, prob = 0.90)

# 3M3. Construct a posterior predictive check for this model and data. This means simulate the distribution of samples, averaging over the posterior uncertainty in p. What is the probability of observing 8 water in 15 tosses?
ppc_15_toss <- rbinom(1e4, size = 15, prob = samples)
mean(ppc_15_toss == 8)

# 3M4. Using the posterior distribution constructed from the new (8/15) data, now calculate the probability of observing 6 water in 9 tosses.
ppc_9_toss <- rbinom(1e4 * 500, size = 9, prob = samples) # very large sample
mean(ppc_9_toss == 6)

mean((dbinom(6, size = 9, prob = samples))) # The above is equivalent to evaluating the likelihood function, marginalizing over the posterior. As the number of PPC sims -> Inf, they give the same result.

# 3M5. Start over at 3M1, but now use a prior that is zero below p = 0.5 and a constant above p = 0.5. This corresponds to prior information that a majority of the Earth’s surface is water. Repeat each problem above and compare the inferences. What difference does the better prior make? If it helps, compare inferences (using both priors) to the true value p = 0.7.
prior2 <- ifelse(p_grid < 0.5, 0, 1)
likelihood <- dbinom(6, size = 9, prob = p_grid)
posterior <- likelihood * prior2
posterior <- posterior / sum(posterior)

samples <- sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)
HPDI(samples, prob = 0.90)

ppc_15_toss <- rbinom(1e4, size = 15, prob = samples)
mean(ppc_15_toss == 8)

ppc_9_toss <- rbinom(1e4 * 500, size = 9, prob = samples) # very large sample
mean(ppc_9_toss == 6)

# 3M6. Suppose you want to estimate the Earth’s proportion of water very precisely. Specifically, you want the 99% percentile interval of the posterior distribution of p to be only 0.05 wide. This means the distance between the upper and lower bound of the interval should be 0.05. How many times will you have to toss the globe to do this?

prior_samp <- sample(p_grid, prob = prior2, size = 1e4, replace = T)

# A function to calculate how wide an interval is
CI_width <- function(toss_sim, n_toss){
  likelihood <- dbinom(toss_sim, size = n_toss, prob = p_grid)
  posterior <- likelihood * prior2
  posterior <- posterior / sum(posterior)
  post_samp <- sample(p_grid, prob = posterior, size = 1e4, replace = T)
  
  return(
    diff(PI(post_samp, prob = 0.99))
    )
}

# A function to give the precision for a given sample size
sample_size_precision <- function(n_toss){
    # simulate from the prior predictive distribution
    toss_sims <- rbinom(500, size = n_toss, prob = prior_samp)
    CIs <- CI_width(toss_sim = toss_sims, n_toss = n_toss)
    interval_width <- max(CIs) # worst case scenario for precision
    
    return(interval_width)
}

# allow us to pass a vector of values as argument
sample_size_precision <- Vectorize(sample_size_precision)
CI_width <- Vectorize(CI_width) 

sample_size <- seq(from = 10, to = 3000, by = 10) # sample sizes to consider

precision <- sample_size_precision(n_toss = sample_size)

plot(precision ~ sample_size, type = "l", lwd = 3)
abline(h = 0.05, lty = "dashed")

needed_sample <- sample_size[precision < 0.05][2] # increase size by one step to accoutn for sim variance
needed_sample

# Validate this approximation
abs(sample_size_precision(needed_sample) - 0.05) < 1e-4

