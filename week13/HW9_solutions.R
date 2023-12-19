library(rethinking)
library(MASS)

# 11E1 If an event has probability 0.35, what are the log-odds of this event?
logit(0.35)

# 11E2 If an event has log-odds 3.2, what is the probability of this event?
inv_logit(3.2)

# Suppose that a coefficient in a logistic regression has value 1.7. What does this imply about the proportional change in odds of the outcome?
exp(1.7) # odds ratio, 5.47 times greater odds

## 11H3
data(eagles)
d <- eagles
d$pirateL <- ifelse( d$P == "L", 1, 0)
d$victimL <- ifelse( d$V == "L", 1, 0)
d$pirateA <- ifelse( d$A == "A", 1, 0)

m_quap <- quap(
  alist(
    y ~ dbinom(n, p),
    logit(p) <- a + bP*pirateL + bV*victimL + bA*pirateA,
    a ~ dnorm(0, 1.5),
    c(bP, bV, bA) ~ dnorm(0, 0.5)
  ),
  data = d,
)

m_ulam <- ulam(
  alist(
    y ~ dbinom(n, p),
    logit(p) <- a + bP*pirateL + bV*victimL + bA*pirateA,
    a ~ dnorm(0, 1.5),
    c(bP, bV, bA) ~ dnorm(0, 0.5)
  ),
  data = d,
  chains = 4,
  cores = 4,
  log_lik = T
)

precis(m_quap)
precis(m_ulam)

# Predicted probability of success for each row of the data
prob_fit <- link(m_ulam)
pred_fit <- sim(m_ulam)

plot(apply(prob_fit, 2, mean), ylim = c(0,1), pch = 16, ylab = "probability of success")
PI_prob <- apply(prob_fit, 2, HPDI)

for (i in 1:nrow(d)) lines(x = rep(i, 2), y = PI_prob[,i])

plot(apply(pred_fit, 2, mean), pch = 16, ylab = "number of successes", ylim = c(0, 30))
PI_pred <- apply(pred_fit, 2, HPDI)

for (i in 1:nrow(d)) lines(x = rep(i, 2), y = PI_pred[,i])

# try interaction model
m_ulam_interaction <- ulam(
  alist(
    y ~ dbinom(n, p),
    logit(p) <- a + bP*pirateL + bV*victimL + bA*pirateA + bPA*pirateL*pirateA,
    a ~ dnorm(0, 1.5),
    c(bP, bV, bA, bPA) ~ dnorm(0, 0.5)
  ),
  data = d,
  chains = 4,
  cores = 4,
  log_lik = T
)

compare(m_ulam, m_ulam_interaction, func = "WAIC")
