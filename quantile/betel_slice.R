# =========================================================
# ETEL(=EL, order 0) likelihood for quantile constraint: E[1(x<=theta)] = p
# =========================================================
# log_likelihood_etel_quantile <- function(theta, x, p) {
#   n <- length(x)
#   k <- sum(x <= theta)  # = n F_n(theta)
#   
#   if (p <= 0 || p >= 1) return(-Inf)
#   
#   if (k == 0) {
#     # ((1-p)/n)^n
#     return(n * (log1p(-p) - log(n)))
#   } else if (k == n) {
#     # (p/n)^n
#     return(n * (log(p) - log(n)))
#   } else {
#     # (p/k)^k * ((1-p)/(n-k))^(n-k)
#     return(k * (log(p) - log(k)) + (n - k) * (log1p(-p) - log(n - k)))
#   }
# }


log_likelihood_etel_quantile <- function(theta, x, p) {
  n <- length(x)
  k <- sum(x <= theta)  # = n F_n(theta)
  
  if (p <= 0 || p >= 1) return(-Inf)
  
  if (k == 0) {
    # ((1-p)/n)^n
    return(log(0))
  } else if (k == n) {
    # (p/n)^n
    return(log(0))
  } else {
    # (p/k)^k * ((1-p)/(n-k))^(n-k)
    return(k * (log(p) - log(k)) + (n - k) * (log1p(-p) - log(n - k)))
  }
}

# =========================================================
# Posterior = ETEL likelihood + Cauchy prior
# =========================================================
log_posterior_etel <- function(theta, x, p, m, tau) {
  log_like <- log_likelihood_etel_quantile(theta, x, p)
  if (!is.finite(log_like)) return(-Inf)
  
  # log_prior <- dcauchy(theta, location = m, scale = tau, log = TRUE)
  log_prior <- dt((theta - m)/tau, df = 2.5, log = TRUE) - log(tau)
  log_like + log_prior
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
