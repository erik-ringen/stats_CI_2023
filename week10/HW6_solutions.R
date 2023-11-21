library(rethinking)

### 7H1
data(Laffer)

Laffer$tax_rate_s <- standardize(Laffer$tax_rate)
Laffer$tax_revenue_s <- standardize(Laffer$tax_revenue)

m_linear <- quap(
  alist(
    tax_revenue_s ~ dnorm(mu, sigma),
    mu <- a + b*tax_rate_s,
    c(a, b) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = Laffer
)

m_quad <- quap(
  alist(
    tax_revenue_s ~ dnorm(mu, sigma),
    mu <- a + b*tax_rate_s + b2*tax_rate_s^2,
    c(a, b, b2) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = Laffer
)

mod_compare <- compare(m_linear, m_quad, func = PSIS, n = 1e4) # pPSIS = penalty term
mod_compare

# Compare the predictions of each model
pred_linear <- link(m_linear, data = Laffer)
pred_quad <- link(m_quad, data = Laffer)

mu_linear <- apply(pred_linear, 2, mean)
mu_quad <- apply(pred_quad, 2, mean)

PI_linear <- apply(pred_linear, 2, PI, prob = 0.9)
PI_quad <- apply(pred_quad, 2, PI, prob = 0.9)

# Plot each model's predictions
par(mfrow = c(1,3))

plot(Laffer$tax_revenue_s ~ Laffer$tax_rate_s, ylab = "tax revenue z-score", xlab = "tax rate z-score")
lines(x = Laffer$tax_rate_s, y = mu_linear)
shade(PI_linear, Laffer$tax_rate_s)
mtext("m_linear")

plot(Laffer$tax_revenue_s ~ Laffer$tax_rate_s, ylab = "tax revenue z-score", xlab = "tax rate z-score")
lines(x = Laffer$tax_rate_s, y = mu_quad)
shade(PI_quad, Laffer$tax_rate_s)
mtext("m_quadratic")

# We could also look at an ensemble of predictions, weighting by the PSIS weights
pred_ensemble <- pred_linear * mod_compare$weight[1] + pred_quad * mod_compare$weight[2] 

mu_ensemble <- apply(pred_ensemble, 2, mean)
PI_ensemble <- apply(pred_ensemble, 2, PI, prob = 0.9)

plot(Laffer$tax_revenue_s ~ Laffer$tax_rate_s, ylab = "tax revenue z-score", xlab = "tax rate z-score")
lines(x = Laffer$tax_rate_s, y = mu_ensemble)
shade(PI_ensemble, Laffer$tax_rate_s)
mtext("Ensemble")

### 7H2
PSIS_quad <- PSIS(m_quad, pointwise = T)
WAIC_quad <- rethinking::WAIC(m_quad, pointwise = T)

plot(PSIS_quad$k ~ WAIC_quad$penalty)
abline(h = 1, lty = "dashed")
which.max(PSIS_quad$penalty) # don't know which country this outlier is

m_linear_robust <- quap(
  alist(
    tax_revenue_s ~ dstudent(2, mu, sigma),
    mu <- a + b*tax_rate_s,
    c(a, b) ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = Laffer
)

m_quad_robust <- quap(
  alist(
    tax_revenue_s ~ dstudent(2, mu, sigma),
    mu <- a + b*tax_rate_s + b2*tax_rate_s^2,
    c(a, b, b2) ~ dnorm(0, 1),
    sigma ~ dexp(1) ),
  data = Laffer
)

mod_compare_robust <- compare(m_linear_robust, m_quad_robust, func = PSIS, n = 1e4) # pPSIS = penalty term
mod_compare_robust

# Compare the predictions of each model
pred_linear_robust <- link(m_linear_robust, data = Laffer)
pred_quad_robust <- link(m_quad_robust, data = Laffer)

mu_linear_robust <- apply(pred_linear_robust, 2, mean)
mu_quad_robust <- apply(pred_quad_robust, 2, mean)

PI_linear_robust <- apply(pred_linear_robust, 2, PI, prob = 0.9)
PI_quad_robust <- apply(pred_quad_robust, 2, PI, prob = 0.9)

# Plot each model's predictions
par(mfrow = c(1,3))

plot(Laffer$tax_revenue_s ~ Laffer$tax_rate_s, ylab = "tax revenue z-score", xlab = "tax rate z-score")
lines(x = Laffer$tax_rate_s, y = mu_linear_robust)
shade(PI_linear_robust, Laffer$tax_rate_s)
mtext("m_linear_robust")

plot(Laffer$tax_revenue_s ~ Laffer$tax_rate_s, ylab = "tax revenue z-score", xlab = "tax rate z-score")
lines(x = Laffer$tax_rate_s, y = mu_quad_robust)
shade(PI_quad_robust, Laffer$tax_rate_s)
mtext("m_quadratic_robust")

# We could also look at an ensemble of predictions, weighting by the PSIS weights
pred_ensemble_robust <- pred_linear_robust * mod_compare_robust$weight[2] + pred_quad_robust * mod_compare_robust$weight[1] 

mu_ensemble_robust <- apply(pred_ensemble_robust, 2, mean)
PI_ensemble_robust <- apply(pred_ensemble_robust, 2, PI, prob = 0.9)

plot(Laffer$tax_revenue_s ~ Laffer$tax_rate_s, ylab = "tax revenue z-score", xlab = "tax rate z-score")
lines(x = Laffer$tax_rate_s, y = mu_ensemble_robust)
shade(PI_ensemble_robust, Laffer$tax_rate_s)
mtext("Ensemble")
