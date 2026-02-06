pkgs <- c("czzg", "RcppXPtrUtils", "betel2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}


Rptr_quant_95 <- cppXPtr('
arma::field<arma::mat> g(const arma::mat theta,
                         const arma::colvec y,
                         const arma::mat dat) {

  int n = y.n_rows;
  
  double q = theta(0,0);
  arma::colvec g = arma::conv_to<arma::colvec>::from(y <= q) - 0.95;

  arma::field<arma::mat> rhof(1,1);
  rhof(0,0) = g;
  return rhof;
}
', depends = "RcppArmadillo")


Rptr_quant_99 <- cppXPtr('
arma::field<arma::mat> g(const arma::mat theta,
                         const arma::colvec y,
                         const arma::mat dat) {

  int n = y.n_rows;
  
  double q = theta(0,0);
  arma::colvec g = arma::conv_to<arma::colvec>::from(y <= q) - 0.99;

  arma::field<arma::mat> rhof(1,1);
  rhof(0,0) = g;
  return rhof;
}
', depends = "RcppArmadillo")


betel_quantile <- function(x, p, psi0, prior_mean, prior_sd, burn, n_samples, seed) {
  n <- length(x)
  k <- length(prior_mean) # number of parameters 
  d <- 1                  # number of moments
  
  Z <- matrix(1, nrow = n, ncol = 1) 
  qZls <- list(Z)
  
  psi0_ <- matrix(prior_mean, ncol=k)
  Psi0_ <- rep(prior_sd^2, k)
  
  if (p == 0.95) {
    rhofunc <- Rptr_quant_95
  } else {
    rhofunc <- Rptr_quant_99
  }
  
  psim <- betel2::bayesetel(
    rhofunc    = rhofunc,
    qZls       = qZls,
    y          = x,
    dat        = matrix(1, nrow = n, ncol = 1), 
    psi0       = psi0,  # start value, sample estimation을 사용 못함.
    lam0       = 0.5*dnorm(1), # start value
    psi0_      = psi0_,
    Psi0_      = Psi0_,
    controlpsi = list(maxiterpsi = 1000, mingrpsi = 1e-6),
    controllam = list(maxiterlam = 100, mingrlam = 1e-5),
    n0         = burn, 
    m          = n_samples - burn,
    printstep  = n_samples, 
    seed       = seed 
  )
  
  psim
}



