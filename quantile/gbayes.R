if (!requireNamespace("nimble", quietly = TRUE)) {
  install.packages("nimble")
}
library(nimble)

estimate_w_joint <- function(x_data,
                             q_levels = c(0.95, 0.99),
                             prior_mean = c(0, 0),  
                             prior_sd = c(10, 10),
                             n_mc = 20000) {
  stopifnot(length(q_levels) == 2)
  
  N <- length(x_data)
  q_emp <- as.numeric(quantile(x_data, q_levels))  
  
  joint_loss <- function(theta, x_subset, qs) {
    if (theta[1] >= theta[2]) return(1e8 + 1e6 * (theta[1]-theta[2])^2)
    l <- 0
    for (j in seq_along(qs)) {
      q <- qs[j]; th <- theta[j]
      l <- l + sum((1 - q) * pmax(th - x_subset, 0) + q * pmax(x_subset - th, 0))
    }
    l
  }
  
  # denominator via LOO-CV
  loo_losses <- numeric(N)
  for (i in 1:N) {
    x_minus_i <- x_data[-i]
    opt <- optim(par = sort(q_emp),
                 fn  = joint_loss,
                 x_subset = x_minus_i,
                 qs = q_levels,
                 method = "Nelder-Mead",
                 control = list(maxit = 5000))
    if (opt$convergence != 0) {
      loo_losses[i] <- NA_real_
      next
    }
    theta_minus_i <- opt$par
    loo_losses[i] <- joint_loss(theta_minus_i, x_data[i], q_levels)
  }
  denom <- mean(loo_losses, na.rm = TRUE)
  
  # numerator: prior information term
  log_prior <- function(theta) {
    if (theta[1] >= theta[2]) return(-Inf)
    sum(dnorm(theta, mean = prior_mean, sd = prior_sd, log = TRUE))
  }
  theta_b <- prior_mean  
  lp_b <- log_prior(theta_b)
  
  samples <- matrix(NA_real_, nrow = n_mc, ncol = 2)
  k <- 0L
  while (k < n_mc) {
    cand <- rnorm(2, mean = prior_mean, sd = prior_sd)
    if (cand[1] < cand[2]) {
      k <- k + 1L
      samples[k, ] <- cand
    }
  }
  lp_s <- apply(samples, 1, log_prior)
  numer <- lp_b - mean(lp_s)
  
  w_hat <- numer / denom
  w_hat
}


estimate_w <- function(x_data,
                       q_level = 0.95,
                       prior_mean = 0,
                       prior_sd = 10,
                       n_mc = 20000) {
  stopifnot(length(q_level) == 1, q_level > 0, q_level < 1)
  
  N <- length(x_data)
  q_emp <- as.numeric(quantile(x_data, q_level)) 
  
  loss_sum <- function(theta, x_subset, q) {
    sum((1 - q) * pmax(theta - x_subset, 0) + q * pmax(x_subset - theta, 0))
  }
  
  # 1) denominator via LOO-CV
  loo_losses <- numeric(N)
  for (i in 1:N) {
    x_minus_i <- x_data[-i]
    opt <- optim(par = q_emp,
                 fn  = loss_sum,
                 x_subset = x_minus_i,
                 q = q_level,
                 method = "Brent",
                 lower = 0.0001, 
                 upper = 1000,
                 control = list(maxit = 5000))
    if (opt$convergence != 0) {
      loo_losses[i] <- NA_real_
      next
    }
    theta_minus_i <- opt$par
    loo_losses[i] <- loss_sum(theta_minus_i, x_data[i], q_level)
  }
  denom <- mean(loo_losses, na.rm = TRUE)
  
  # 2) numerator: prior information term 
  log_prior <- function(theta) dnorm(theta, mean = prior_mean, sd = prior_sd, log = TRUE)
  
  theta_b <- prior_mean                  
  lp_b <- log_prior(theta_b)
  
  prior_samples <- rnorm(n_mc, mean = prior_mean, sd = prior_sd)
  lp_s <- dnorm(prior_samples, mean = prior_mean, sd = prior_sd, log = TRUE)
  
  numer <- lp_b - mean(lp_s)
  
  w_hat <- numer / denom
  w_hat
}




JointQuantileCode <- nimbleCode({
  theta95 ~ dnorm(prior_mean[1], sd = prior_sd[1])
  theta99 ~ dnorm(prior_mean[2], sd = prior_sd[2])

  constraint ~ dconstraint(theta95 < theta99)
  
  for(i in 1:N) {
    l95[i] <- 0.05 * step(theta95 - x[i]) * (theta95 - x[i]) + 0.95 * step(x[i] - theta95) * (x[i] - theta95)
    l99[i] <- 0.01 * step(theta99 - x[i]) * (theta99 - x[i]) + 0.99 * step(x[i] - theta99) * (x[i] - theta99)
    
    totalLoss[i] <- w * (l95[i] + l99[i])
    zeros[i] ~ dpois(totalLoss[i])
  }
})


QuantileCode <- nimbleCode({
  theta ~ dnorm(prior_mean, sd = prior_sd)

  for(i in 1:N) {
    loss[i] <- (1 - p) * step(theta - x[i]) * (theta - x[i]) + p * step(x[i] - theta) * (x[i] - theta)
    
    totalLoss[i] <- w * loss[i]
    zeros[i] ~ dpois(totalLoss[i])
  }
})


gbayes_quantile <- function(x, p, prior_mean, prior_sd, burn, n_samples, nb_code = QuantileCode) {
  w_est <- estimate_w(x_data = x, q_level = p, n_mc = 5000, 
                      prior_mean = prior_mean, prior_sd = prior_sd)
  
  constants <- list(N = length(x), w = w_est, p = p, prior_mean = prior_mean, prior_sd = prior_sd)
  data <- list(x = x, zeros = rep(0, length(x))) 
  inits <- list(theta = quantile(x, probs = p, names = FALSE))
  
  model <- nimbleModel(nb_code, constants, data, inits)
  cModel <- compileNimble(model)
  mcmcConf <- configureMCMC(model)
  mcmc <- buildMCMC(mcmcConf)
  cMcmc <- compileNimble(mcmc, project = model)
  psamp <- runMCMC(cMcmc, niter = n_samples, nburnin = burn)
  nimble::clearCompiled(model)
  gc()
  psamp
}