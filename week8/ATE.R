library(rethinking)
library(ggplot2)

set.seed(2)
N <- 3000
Z <- sample(0:1, size = N, replace = T)
X <- rnorm(N)
Y <- rnorm(N, ifelse(Z == 0, 1*X, -1*X))

d <- data.frame(X = X, Z = as.factor(Z), Y = Y)

ggplot(d, aes(x = X, y = Y, color = Z)) + geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", se = F) + 
  theme_minimal(base_size = 18)

ggplot(d, aes(x = X, y = Y)) + geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", se = F) + 
  theme_minimal(base_size = 18)

# Fit a model ignoring effect heterogeneity
m1 <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + b*X,
    c(a, b) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = d
)

summary(m1)

# Fit a model that stratifies by Z
d_list <- list(
  Y_0 = Y[Z == 0],
  Y_1 = Y[Z == 1],
  X_0 = X[Z == 0],
  X_1 = X[Z == 1]
)

m2 <- quap(
  alist(
    Y_0 ~ dnorm(mu_0, sigma),
    Y_1 ~ dnorm(mu_1, sigma),
    
    mu_0 <- a_0 + b_0*X_0,
    mu_1 <- a_1 + b_1*X_1,
    c(a_0, a_1, b_0, b_1) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = d_list
)

# Calculate the average treatment effect from each model
post_m1 <- extract.samples(m1)
post_m2 <- extract.samples(m2)

# What if we had intervened so that X = 0 and X = 1?
P_Y0_m1 <- sapply(rep(0, N), function(x) rnorm(1e4, post_m1$a + post_m1$b*x, post_m1$sigma)) 
P_Y1_m1 <- sapply(rep(1, N), function(x) rnorm(1e4, post_m1$a + post_m1$b*x, post_m1$sigma)) 

ATE_m1 <- apply( P_Y1_m1 - P_Y0_m1, 2, mean )
dens(ATE_m1)

# For the model that stratifies by Z, we must also stratify our predictions
P_Y0_m2_Z0 <- sapply(rep(0, N), function(x) rnorm(1e4, post_m2$a_0 + post_m2$b_0*x, post_m2$sigma)) 
P_Y1_m2_Z0 <- sapply(rep(1, N), function(x) rnorm(1e4, post_m2$a_0 + post_m2$b_0*x, post_m2$sigma)) 

P_Y0_m2_Z1 <- sapply(rep(0, N), function(x) rnorm(1e4, post_m2$a_1 + post_m2$b_1*x, post_m2$sigma)) 
P_Y1_m2_Z1 <- sapply(rep(1, N), function(x) rnorm(1e4, post_m2$a_1 + post_m2$b_1*x, post_m2$sigma)) 

# Different treatment effects depending on Z
ATE_m2_Z0 <- apply( P_Y1_m2_Z0 - P_Y0_m2_Z0, 2, mean )
ATE_m2_Z1 <- apply( P_Y1_m2_Z1 - P_Y0_m2_Z1, 2, mean )

# Finally, for the population average treatment effect, we weight by Pr(Z)
Pr_Z1 <- mean(Z == 1)

ATE_m2 <- ATE_m2_Z0 * (1 - Pr_Z1) + ATE_m2_Z1 * (Pr_Z1)

# Now compare ATE from both models
dens(ATE_m1, ylim = c(0,50), col = "skyblue")
dens(ATE_m2, add = T, col = "red")
