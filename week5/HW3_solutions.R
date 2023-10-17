library(rethinking)

#### 4H1
data(Howell1)

d <- Howell1
d2 <- d[d$age >= 18, ]

m <- quap(
    alist(
      height ~ dnorm(mu,sigma),
      mu <- a + b*weight,
      a ~ dnorm(178, 20),
      b ~ dnorm(0,10),
      sigma ~ dunif(0,50)),
    data = d2
)

# extract posterior samps
post <- extract.samples(m)

new_weights <- c(46.95, 43.72, 64.78, 32.59, 54.63)

# Simulate heights for these individuals
new_heights <- matrix(NA, nrow = 1e4, ncol = length(new_weights))

for (j in 1:length(new_weights)) {
  new_heights[,j] = rnorm(1e4, post$a + post$b*new_weights[j], post$sigma)
}

head(new_heights)

# Expected heights for each individual
apply(new_heights, 2, mean)

# 89% HPDI for each individual
apply(new_heights, 2, HPDI, prob = 0.89)

#### 4H2
d_kids <- d[d$age < 18, ]

m_kids <- quap(
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b*weight,
    a ~ dnorm(108, 20),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0,50)),
  data = d_kids
)

post <- extract.samples(m_kids)

# For every 10 units increase in weight, how much taller does the model predict a child will get?
mean( post$b * 10 )
HPDI( post$b * 10 )

# Plot data vs model preds
weight_seq <- seq(from = min(d_kids$weight), to = max(d_kids$weight), length.out = 30)

# expected values
mu_height <- sapply(weight_seq, function(x) post$a + post$b*x)
# predicted values
pred_height <- sapply(weight_seq, function(x) rnorm(1e4, post$a + post$b*x, post$sigma))

# plot data
with(d_kids, plot(height ~ weight, pch = 16, col = col.alpha("cornflowerblue", 0.4)))

lines(weight_seq, apply(mu_height, 2, mean), lwd = 2)
shade(apply(mu_height, 2, HPDI), weight_seq)
shade(apply(pred_height, 2, HPDI), weight_seq, col = col.alpha("black", 0), border = "black", lty = "dashed")

# Check residuals
model_epred <- sapply(d_kids$weight, function(x) post$a + post$b*x)

residual <- d_kids$height - apply(model_epred, 2, mean)

# check for normality of residuals
dens(residual, norm.comp = T)
qqnorm(residual, pch = 16, col = col.alpha("black", 0.4))
qqline(residual, col = "steelblue", lwd = 2)

# Generic diagnostics can tell us the model has a problem, but not where the problem is coming from. More useful to examine deviations from the linear model expectation
plot(residual ~ d_kids$weight, pch = 16, col = col.alpha("black", 0.4))
abline(h = 0, lty = "dashed")

#### 4H3
Howell1$log_weight <- log(Howell1$weight)

m_log <- quap(
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b*log_weight,
    a ~ dnorm(138, 20),
    b ~ dnorm(0, 10),
    sigma ~ dunif(0,50)),
  data = Howell1
)

post <- extract.samples(m_log)

log_weight_seq <- seq(from = min(Howell1$log_weight), to = max(Howell1$log_weight), length.out = 30)

# expected values
mu_height <- sapply(log_weight_seq, function(x) post$a + post$b*x)
# predicted values
pred_height <- sapply(log_weight_seq, function(x) rnorm(1e4, post$a + post$b*x, post$sigma))

# On log scale (log-linear)
with(Howell1, plot(height ~ log(weight), pch = 16, col = col.alpha("cornflowerblue")))

lines(log_weight_seq, apply(mu_height, 2, mean))
shade(apply(mu_height, 2, HPDI, prob = 0.97), log_weight_seq)
shade(apply(pred_height, 2, HPDI, prob = 0.97), log_weight_seq, col = col.alpha("black", 0), border = "black", lty = "dashed")

# On original measurement scale
with(Howell1, plot(height ~ weight, pch = 16, col = col.alpha("cornflowerblue")))

lines(exp(log_weight_seq), apply(mu_height, 2, mean))
shade(apply(mu_height, 2, HPDI, prob = 0.97), exp(log_weight_seq))
shade(apply(pred_height, 2, HPDI, prob = 0.97), exp(log_weight_seq), col = col.alpha("black", 0), border = "black", lty = "dashed")

#### 4H4
d <- Howell1

d$weight_s <- (d$weight - mean(d$weight)) / sd(d$weight)
d$weight_s2 <- d$weight_s^2

# Original model from text
m4.5 <- quap(
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + b1*weight_s + b2*weight_s2,
    a ~ dnorm(178,20),
    b1 ~ dlnorm(0,1),
    b2 ~ dnorm(0,1),
    sigma ~ dunif(0,50)
  ) , data = d)

# prior predictions for height
a <- rnorm(1e4, 178, 20)
b1 <- rlnorm(1e4, 0, 1)
b2 <- rnorm(1e4, 0, 1)
sigma <- runif(1e4, 0, 50)

prior_pred_height <- sapply(d$weight_s, function(x) rnorm(1e4, a + b1*x + b2*x^2, sigma))

dens(prior_pred_height, ylim = c(0, 0.1)) # prior predictive distribution
for (i in 1:20) dens(prior_pred_height[i,], col = col.alpha("black"), add = T) # many possible realizations of the data
dens(d$height, col = "darkred", add = T, lwd = 2)

# prior distribution of mean height
dens(apply(prior_pred_height, 1, mean))
abline(v = mean(d$height), col = "darkred")

# prior distribution of min height
dens(apply(prior_pred_height, 1, min))
abline(v = min(d$height), col = "darkred")

# prior distribution of max height
dens(apply(prior_pred_height, 1, max))
abline(v = max(d$height), col = "darkred")

# Prior distribution of functions
curve(a[1] + b1[1]*x + b2[1]*x^2, from = min(d$weight_s), to = max(d$weight_s), ylim = c(min(d$height)*0.5, max(d$height)*1.5), col = col.alpha("black"))

for (i in 1:200) {
  curve(a[i] + b1[i]*x + b2[i]*x^2, from = min(d$weight_s), to = max(d$weight_s), col = col.alpha("black"), add = T)
}

points(d$height ~ d$weight_s, col = "darkred")

# trying to better with diff priors
a <- rnorm(1e4, mean(d$height), 20) # center on sample mean
b1 <- abs(rnorm(1e4, 0, 10)) # allow for larger linear effect
b2 <- rnorm(1e4, 0, 5) # allow for larger quadratic effect
sigma <- runif(1e4, 0, 10) # reduce std deviation for less variance in predictions


prior_pred_height <- sapply(d$weight_s, function(x) rnorm(1e4, a + b1*x + b2*x^2, sigma))

dens(prior_pred_height, ylim = c(0, 0.1)) # prior predictive distribution
for (i in 1:20) dens(prior_pred_height[i,], col = col.alpha("black"), add = T) # many possible realizations of the data
dens(d$height, col = "darkred", add = T, lwd = 2)

# prior distribution of mean height
dens(apply(prior_pred_height, 1, mean))
abline(v = mean(d$height), col = "darkred")

# prior distribution of min height
dens(apply(prior_pred_height, 1, min))
abline(v = min(d$height), col = "darkred")

# prior distribution of max height
dens(apply(prior_pred_height, 1, max))
abline(v = max(d$height), col = "darkred")

# Prior distribution of functions
curve(a[1] + b1[1]*x + b2[1]*x^2, from = min(d$weight_s), to = max(d$weight_s), ylim = c(min(d$height)*0.5, max(d$height)*1.5), col = col.alpha("black"))

for (i in 1:200) {
  curve(a[i] + b1[i]*x + b2[i]*x^2, from = min(d$weight_s), to = max(d$weight_s), col = col.alpha("black"), add = T)
}

points(d$height ~ d$weight_s, col = "darkred")












