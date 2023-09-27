# Analytically update the posterior distribution
globe_toss <- function(prior_alpha = 1, prior_beta = 1, W, L) {
  post_alpha <- prior_alpha + W
  post_beta <- prior_beta + L
  
  # plot density
  curve(dbeta(x, post_alpha, post_beta), from = 0, to = 1, ylab = "Probability density", xlab = "p (proportion of water)")
  mtext(paste("W =", W, " L =", L))
}

# update the values of W and L according to the observed data
globe_toss(W = 9, L = 2)

