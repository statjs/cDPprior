if (!requireNamespace("nimble", quietly = TRUE)) {
  install.packages("nimble")
}
library(nimble)


estimate_w_reg <- function(y, X, prior_mean = 0, prior_sd = 10, n_mc = 5000) {
  N <- length(y)
  
  # Squared Error Loss
  loss_reg <- function(beta, y_sub, X_sub) {
    sum((y_sub - X_sub %*% beta)^2)  # e = y - X * beta
  }
  
  # LOO-CV
  loo_losses <- numeric(N)
  beta_init <- solve(t(X) %*% X) %*% t(X) %*% y 
  
  for (i in 1:N) {
    opt <- optim(par = beta_init, fn = loss_reg, 
                 y_sub = y[-i], X_sub = X[-i, , drop=FALSE],
                 method = "BFGS") 
    
    if (opt$convergence == 0) {
      loo_losses[i] <- loss_reg(opt$par, y[i], X[i, , drop=FALSE])
    } else {
      loo_losses[i] <- NA_real_
    }
  }
  denom <- mean(loo_losses, na.rm = TRUE)
  
  # prior
  log_prior <- function(beta) sum(dnorm(beta, mean = prior_mean, sd = prior_sd, log = TRUE))
  numer <- log_prior(rep(prior_mean, ncol(X))) - 
    mean(replicate(n_mc, log_prior(rnorm(ncol(X), prior_mean, prior_sd))))
  
  return(abs(numer / denom))
}


# Regression - Nimble
MeanRegCode <- nimbleCode({
  for(k in 1:K) {
    beta[k] ~ dnorm(prior_mean, sd = prior_sd)
  }
  
  for(i in 1:N) {
    mu[i] <- inprod(X[i, 1:K], beta[1:K])
    loss[i] <- (y[i] - mu[i])^2
    
    totalLoss[i] <- w * loss[i]
    zeros[i] ~ dpois(totalLoss[i])
  }
})


gbayes_regression <- function(y, X, prior_mean = 0, prior_sd = 10, 
                              burn = 5000, n_samples = 10000) {
  
  w_est <- estimate_w_reg(y, X, prior_mean = prior_mean, prior_sd = prior_sd)
  
  constants <- list(N = length(y), K = ncol(X), w = w_est, 
                    prior_mean = prior_mean, prior_sd = prior_sd)
  data <- list(y = as.numeric(y), X = X, zeros = rep(0, length(y)))
  inits <- list(beta = as.vector(solve(t(X) %*% X) %*% t(X) %*% y))
  
  model <- nimbleModel(MeanRegCode, constants, data, inits)
  cModel <- compileNimble(model)
  mcmcConf <- configureMCMC(model)
  mcmc <- buildMCMC(mcmcConf)
  cMcmc <- compileNimble(mcmc, project = model)
  psamp <- runMCMC(cMcmc, niter = n_samples, nburnin = burn)
  nimble::clearCompiled(model)
  gc()
  psamp
}
