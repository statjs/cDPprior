pkgs <- c("czzg", "RcppXPtrUtils", "betel2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

Rmatptr = cppXPtr('
arma::field<arma::mat> g(const arma::mat theta,
                         const arma::colvec y,
                         const arma::mat dat) {
  
  arma::colvec e = y - dat * theta;
  
  arma::field<arma::mat> rhof(1, 1);
  rhof(0, 0) = e;               

  return rhof;
}
', depends = "RcppArmadillo")


betel_regression <- function(y, X, psi0, prior_mean, prior_sd, burn, n_samples, seed) {
  n <- length(y)
  k <- ncol(X) 
  
  Z <- cbind(1, X)  
  qZls <- vector("list", 1) 
  qZls[[1]] <- Z
  d <- ncol(qZls[[1]])

  psi0_ <- matrix(prior_mean, nrow = k, ncol = 1)
  Psi0_ <- rep(prior_sd, k)
  
  set.seed(1)
  psim <- betel2::bayesetel(
    rhofunc    = Rmatptr,
    qZls       = qZls,
    y          = y,
    dat        = cbind(X),
    psi0       = psi0,
    lam0       = 0.5 * rnorm(d),
    psi0_      = psi0_,
    Psi0_      = Psi0_,
    controlpsi = list(maxiterpsi = 1000, mingrpsi = 1e-7),
    controllam = list(maxiterlam = 100, mingrlam = 1e-6),
    n0         = burn, 
    m          = n_samples - burn, 
    printstep  = n_samples,  
    seed       = seed
  )
  
  psim
}
