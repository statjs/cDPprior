# ==============================================================================
# Jaeyong Lee, Kwangmin Lee, Jaegui Lee, and Seongil Jo (2025).
# Conditional Dirichlet Processes and Functional Condition Models.
# arXiv:2506.15932
# ==============================================================================
## Moments Estimation
## Figure S1 & Table S1

## -- packages
library(tidyverse)
library(boot)
library(e1071)  # skewness 

## -- somites data
somites <- readr::read_csv("somites_data.csv")
somite_data <- rep(somites$somites, somites$frequencies)


# ==============================================================================
# Figure S1
df <- data.frame(somite = somite_data)

ggplot(df, aes(x = somite)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(x = "Number of Somites", y = "Frequency") +
  theme_minimal()


# ==============================================================================
# Table S1

set.seed(1)
(boot_mean <- boot(data = somite_data, statistic = function(x, i) mean(x[i]), R = 1000))
boot.ci(boot_mean, type = "perc")

(boot_sd <- boot(data = somite_data, statistic = function(x, i) sd(x[i]), R = 1000))
boot.ci(boot_sd, type = "perc")

(boot_skew <- boot(data = somite_data, statistic = function(x, i) skewness(x[i]), R = 1000))
boot.ci(boot_skew, type = "perc")

(boot_kurt <- boot(data = somite_data, statistic = function(x, i) kurtosis(x[i]), R = 1000))
boot.ci(boot_kurt, type = "perc")


# ==============================================================================
# cDP prior
fit_cDP <- cDPprior::cdp_4moments(x = somite_data, n_samples = 10000, sep = TRUE, total_mass = 1, N = 100, seed = 131)
summary(fit_cDP, burn = 5000)
