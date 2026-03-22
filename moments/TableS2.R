# ==============================================================================
# Jaeyong Lee, Kwangmin Lee, Jaegui Lee, and Seongil Jo (2025).
# Conditional Dirichlet Processes and Functional Condition Models.
# arXiv:2506.15932
# ==============================================================================
## Moments Estimation under IQR(F) = c
## Table S2

## -- somites data
somites <- readr::read_csv("somites_data.csv")
somite_data <- rep(somites$somites, somites$frequencies)
iqr <- IQR(somite_data)

# ==============================================================================
# cDP prior
fit_cDP <- cDPprior::cdp_4moments(x = somite_data, 
                                  iqr = 9,
                                  n_samples = 10000, 
                                  sep = FALSE, 
                                  total_mass = 1, N = 100, 
                                  seed = 137)
summary(fit_cDP, burn = 5000)
