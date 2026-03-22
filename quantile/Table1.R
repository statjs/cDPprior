rm(list=ls())
# ==============================================================================
# Jaeyong Lee, Kwangmin Lee, Jaegui Lee, and Seongil Jo (2025).
# Conditional Dirichlet Processes and Functional Condition Models.
# arXiv:2506.15932
# ==============================================================================
## Quantile Estimation
## Table 1

## -- packages
pkgs <- c("xtable", "tidyverse")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}


## -- functions
source("SampQuant.R")
source("gbayes.R")
source("betel_slice.R")

mcmc_summary <- function(x) {
  mcmc_obj <- as.mcmc(x)
  psout <- summary(mcmc_obj)
  ess_val <- coda::effectiveSize(mcmc_obj)
  
  c(psout$statistics[1], 
    psout$quantiles['2.5%'],
    psout$quantiles['97.5%'],
    ESS = as.numeric(ess_val))
}

## -- hurricane data
hurricane <- readr::read_csv("HourlyMaxWinds18512006.csv")
x <- hurricane %>%
  filter(Yr == 2006) %>% select(Wmax) %>% 
  unlist() %>% as.numeric()
n <- length(x)

## -- mcmc parameters
burn <- 5000
n_samples <- 10000 # includes the burn-in period

### ============================================================================
### p = 0.95 # 95th quantile
p <- 0.95

# sample quantile
(qhat_sq <- sample_quantile(x, p = p))

# ==============================================================================
# conditional DP (Lee, Lee, Lee, and Jo, 2026)
# priors
prior_sd <- 100
prior_mean <- round(quantile(x, probs = p, names = FALSE))

# cDP prior & cauchy
fout_cDP <- cDPprior::cdp_quantile(x, p = p, 
                                   base_dist = 1, 
                                   prior_par = list(m = prior_mean),
                                   n_samples = n_samples, seed = 0950)
(qhat_cdp <- mcmc_summary(fout_cDP$samples[(burn+1):n_samples]))


# cDP prior & normal
fout_cDP_nrm <- cDPprior::cdp_quantile(x, p = p, 
                                       base_dist = 2, 
                                       prior_par = list(m = prior_mean),
                                       n_samples = n_samples, seed = 0951)
(qhat_cdp_nrm <- mcmc_summary(fout_cDP_nrm$samples[(burn+1):n_samples]))


# ==============================================================================
# generalized Bayes (Bissiri, Holmes, and Walker, 2016)
set.seed(0952)
fout_gb <- gbayes_quantile(x, p = p, 
                           prior_mean = prior_mean, prior_sd = prior_sd, 
                           burn = burn, n_samples = n_samples, 
                           nb_code = QuantileCode)
(qhat_gb <- mcmc_summary(fout_gb))


# ==============================================================================
# BETEL 
set.seed(0954)
fout_betel <- slice_sampler(log_posterior = log_posterior_etel,
                            theta0 = quantile(x, p),
                            w = 0.5,
                            n_samples = n_samples,
                            max_step_out = 10,
                            x = x, p = p, 
                            m = 0, 
                            tau = 5)
(qhat_betel <-   mcmc_summary(fout_betel[(burn+1):n_samples]))


# ==============================================================================
# Table
(tab1.1 <- rbind(SampQuant = c(qhat_sq, NA), 
                 gBayes = qhat_gb, 
                 ETEL = qhat_betel,
                 cDP = qhat_cdp,
                 cDP_Nrm = qhat_cdp_nrm))


### ============================================================================
### p = 0.95 # 95th quantile
p <- 0.99

# sample quantile
(qhat_sq <- sample_quantile(x, p = p))

# ==============================================================================
# conditional DP (Lee, Lee, Lee, and Jo, 2026)
# priors
prior_sd <- 100
prior_mean <- round(quantile(x, probs = p, names = FALSE))

# cDP prior & cauchy
fout_cDP <- cDPprior::cdp_quantile(x, p = p, 
                                   base_dist = 1, 
                                   prior_par = list(m = prior_mean),
                                   n_samples = n_samples, seed = 0990)
(qhat_cdp <- mcmc_summary(fout_cDP$samples[(burn+1):n_samples]))


# cDP prior & normal
fout_cDP_nrm <- cDPprior::cdp_quantile(x, p = p, 
                                       base_dist = 2, 
                                       prior_par = list(m = prior_mean),
                                       n_samples = n_samples, seed = 0991)
(qhat_cdp_nrm <- mcmc_summary(fout_cDP_nrm$samples[(burn+1):n_samples]))


# ==============================================================================
# generalized Bayes (Bissiri, Holmes, and Walker, 2016)
set.seed(0992)
fout_gb <- gbayes_quantile(x, p = p, 
                           prior_mean = prior_mean, prior_sd = prior_sd, 
                           burn = burn, n_samples = n_samples, 
                           nb_code = QuantileCode)
(qhat_gb <- mcmc_summary(fout_gb))



# ==============================================================================
# BETEL 
set.seed(0994)
fout_betel <- slice_sampler(log_posterior = log_posterior_etel,
                            theta0 = quantile(x, p),
                            w = 0.5,
                            n_samples = n_samples,
                            max_step_out = 10,
                            x = x, p = p, 
                            m = 0, 
                            tau = 5)
(qhat_betel <-   mcmc_summary(fout_betel[(burn+1):n_samples]))


# ==============================================================================
# Table
(tab1.2 <- rbind(SampQuant = c(qhat_sq, NA), 
                 gBayes = qhat_gb, 
                 ETEL = qhat_betel,
                 cDP = qhat_cdp,
                 cDP_Nrm = qhat_cdp_nrm))


# ==============================================================================
# Table
(tab1.1 <- as.data.frame(tab1.1))
(tab1.2 <- as.data.frame(tab1.2))

tab1.1$Quantile <- "0.95"
tab1.2$Quantile <- "0.99"

tab1.1$Method <- rownames(tab1.1)
tab1.2$Method <- rownames(tab1.2)

rownames(tab1.1) <- NULL
rownames(tab1.2) <- NULL

tab_all <- rbind(tab1.1, tab1.2)
colnames(tab_all) <- c('Est', '2.5%', '97.5%', 'ESS', 'Quantile', 'Method')
tab_all <- tab_all[, c("Quantile", "Method", "Est", "2.5%", "97.5%", "ESS")]

latex_tab <- xtable(
  tab_all,
  caption = "Comparison of 0.95 and 0.99 Quantile Estimates",
  label   = "tab:quantile_comparison",
  digits  = c(1, 1, 1, 1, 0, 0, 0)
)

print(
  latex_tab,
  include.rownames = FALSE,
  sanitize.text.function = identity
)
