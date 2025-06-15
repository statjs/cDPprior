# cauchy prior
log_posterior_cauchy <- function(theta, x, p, mu0, sigma0, m, tau) {
  n_theta <- sum(x <= theta)
  n_total <- length(x)
  # alpha_s <- pnorm(theta, mean = mu0, sd = sigma0)
  alpha_s <- pcauchy(theta, location = mu0, scale = sigma0)
  alpha_sc <- 1 - alpha_s
  
  # rising factorial 
  log_rising <- function(a, k) {
    if (k == 0) return(0)
    if (a <= 0) return(-Inf)
    sum(log(a + 0:(k - 1)))
  }
  
  if (alpha_s <= 0 || alpha_s >= 1) return(-Inf)
  
  log_likelihood <- 
    n_theta * log(p) + (n_total - n_theta) * log(1-p) -
    log_rising(alpha_s, n_theta) -
    log_rising(alpha_sc, n_total - n_theta)
  
  log_prior <- dcauchy(theta, location = m, scale = tau, log = TRUE)
  
  return(log_likelihood + log_prior)
}


slice_sampler <- function(log_posterior, theta0, w = 1, n_samples = 5000, max_step_out = 10, ...) {
  samples <- numeric(n_samples)
  samples[1] <- theta0
  
  for (i in 2:n_samples) {
    theta_curr <- samples[i - 1]
    log_post_curr <- log_posterior(theta_curr, ...)
    log_y <- log_post_curr - rexp(1)
    
    u <- runif(1)
    L <- theta_curr - w * u
    R <- L + w
    
    # limit stepping out
    k <- 0
    while (log_posterior(L, ...) > log_y && k < max_step_out) {
      L <- L - w
      k <- k + 1
    }
    k <- 0
    while (log_posterior(R, ...) > log_y && k < max_step_out) {
      R <- R + w
      k <- k + 1
    }
    
    # shrinkage
    repeat {
      theta_new <- runif(1, L, R)
      log_post_new <- log_posterior(theta_new, ...)
      if (log_post_new >= log_y) {
        samples[i] <- theta_new
        break
      } else {
        if (theta_new < theta_curr) {
          L <- theta_new
        } else {
          R <- theta_new
        }
      }
    }
  }
  return(samples)
}
