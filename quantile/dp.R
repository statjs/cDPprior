## =========================================================
## 1) Stick-breaking DP draw (truncated at K)
## =========================================================
rdp <- function(total_mass, base_pr_measure, K = 1000, ...) {
  stopifnot(K >= 2)
  
  betas  <- rbeta(K - 1, 1, total_mass)
  remain <- cumprod(c(1, 1 - betas))          # length K
  
  w <- numeric(K)
  w[1:(K - 1)] <- betas * remain[1:(K - 1)]
  w[K] <- remain[K]                          # exact remainder
  
  atoms <- base_pr_measure(K, ...)
  list(weights = w, atoms = atoms)
}

## =========================================================
## 2) Sample atoms from posterior base G_n
##    G_n = (alpha G0 + sum delta_x) / (alpha + n)
## =========================================================
rG_n <- function(m, x, alpha, rG0, ...) {
  n <- length(x)
  u <- runif(m)
  out <- numeric(m)
  
  take0 <- (u < alpha / (alpha + n))    # from G0
  k0 <- sum(take0)
  
  if (k0 > 0) out[take0] <- rG0(k0, ...)
  if (k0 < m) out[!take0] <- sample(x, size = m - k0, replace = TRUE)
  
  out
}

## =========================================================
## 3) Weighted quantiles (multiple probs in one sort)
## =========================================================
wquantile <- function(atoms, weights, p) {
  atoms <- as.numeric(atoms)
  weights <- as.numeric(weights)
  
  ok <- is.finite(atoms) & is.finite(weights)
  atoms <- atoms[ok]; weights <- weights[ok]
  
  o <- order(atoms)
  atoms <- atoms[o]
  weights <- weights[o]
  
  # safety
  weights <- pmax(weights, 0)
  sw <- sum(weights)
  if (!is.finite(sw) || sw <= 0) return(NA_real_)
  weights <- weights / sw
  
  cw <- cumsum(weights)
  cw[length(cw)] <- 1  # 핵심: 마지막을 1로 고정
  
  idx <- which(cw >= p)[1]
  if (is.na(idx)) idx <- length(atoms)
  atoms[idx]
}

wquantiles <- function(atoms, weights, p = c(0.25, 0.5, 0.75)) {
  sapply(p, function(p) wquantile(atoms, weights, p))
}

## =========================================================
## 4) DP posterior draws for median + IQR 
## =========================================================
dp_post_median_iqr <- function(x,
                               n_samples = 10000,
                               alpha = 1,
                               K = 100,
                               rG0 = function(m, location, scale) rcauchy(m, location, scale),
                               m0 = median(x),
                               s0 = IQR(x) / 2) {
  n <- length(x)
  alpha_post <- alpha + n
  
  out <- data.frame(q25 = numeric(n_samples), 
                    q50 = numeric(n_samples),
                    q75 = numeric(n_samples), 
                    iqr = numeric(n_samples))
  
  for (s in 1:n_samples) {
    dr <- rdp(alpha_post, rG_n, K = K, x = x, alpha = alpha, rG0 = rG0, location = m0, scale = s0)
    
    w <- dr$weights
    a <- dr$atoms
    
    qs <- wquantiles(a, w, p = c(0.25, 0.5, 0.75))
    out$q25[s]  <- qs[1]
    out$q50[s]  <- qs[2]
    out$q75[s]  <- qs[3]
    out$iqr[s]  <- out$q75[s] - out$q25[s]
  }
  
  out
}