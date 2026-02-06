sample_quantile <- function(x, p, alpha = 0.05) {
  n <- length(x)
  qhat <- as.numeric(stats::quantile(x, probs = p, names = FALSE))
  
  # density estimate at qhat
  dens <- stats::density(x)
  fhat <- stats::approx(dens$x, dens$y, xout = qhat, rule = 2)$y
  fhat <- max(fhat, 1e-12)
  
  se <- sqrt(p * (1 - p) / (n * fhat^2))
  z  <- stats::qnorm(1 - alpha / 2)
  
  c(est = qhat, lower = qhat - z * se, upper = qhat + z * se)
}