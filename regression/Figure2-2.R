rm(list=ls())
# ==============================================================================
# Jaeyong Lee, Kwangmin Lee, Jaegui Lee, and Seongil Jo (2025).
# Conditional Dirichlet Processes and Functional Condition Models.
# arXiv:2506.15932
# ==============================================================================
## Regression model
## Figure 2-2

## -- packages
pkgs <- c("xtable", "tidyverse")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}


## -- functions
source("mle.R")
source("gbayes.R")
source("betel.R")

mcmc_summary <- function(x) {
  mcmc_obj <- as.mcmc(x)
  psout <- summary(mcmc_obj)
  ess_val <- coda::effectiveSize(mcmc_obj)
  
  cbind(est = psout$statistics[,1], 
        `2.5%` = psout$quantiles[,'2.5%'],
        `97.5%` = psout$quantiles[,'97.5%'],
        ESS = as.numeric(ess_val))
}


## -- nhanes data
load("nhanes.RData")

y <- nhanes_clean %>% select(LBXGLU) %>% unlist() %>% as.numeric()
yc <- scale(y, scale = FALSE)
n <- length(y)

X <- nhanes_clean %>% select(-c(LBXGLU, SEQN, BMXHT, RIAGENDR)) %>% as.matrix()
Xc <- scale(X, scale = FALSE)
cnames <- colnames(X)
p = dim(Xc)[2]

Zc <- cbind(yc, Xc)

## -- mcmc parameters
burn <- 5000
n_samples <- 10000 # includes the burn-in period

## -- priors
prior_sd <- 10
prior_mean <- 0

## =============================================================================
## mle regression
fout_mle <- mle_regression(y = yc, X = Xc)
(bhat_mle <- cbind(fout_mle, ESS = NA))

## =============================================================================
## conditional Dirichlet process regression
fout_cdp <- cDPprior::cdp_regress(y = yc, X = Xc, n_samples = n_samples, 
                                  base_dist = 1, 
                                  base_par = diag(apply(Zc, 2, IQR)/2),
                                  prior_par = list(m = rep(prior_mean, p), 
                                                   tau = rep(prior_sd, p)),
                                  seed = 771)
(bhat_cdp <- mcmc_summary(fout_cdp$samples[(burn+1):n_samples,]))


## conditional Dirichlet process regression
fout_cdp_nrm <- cDPprior::cdp_regress(y = yc, X = Xc, n_samples = n_samples, 
                                      base_dist = 2, 
                                      base_par = diag(apply(Zc, 2, IQR)/2),
                                      prior_par = list(m = rep(prior_mean, p), 
                                                       tau = rep(prior_sd, p)),
                                      seed = 774)
(bhat_cdp_nrm <- mcmc_summary(fout_cdp_nrm$samples[(burn+1):n_samples,]))

## =============================================================================
## gbayes regression
set.seed(772)
fout_gb <- gbayes_regression(y = yc, X = Xc, 
                             prior_mean = prior_mean, prior_sd = prior_sd, 
                             burn = burn, n_samples = n_samples)
(bhat_gb <- mcmc_summary(fout_gb))


## =============================================================================
## betel regression
fout_bt <- tryCatch({
  betel_regression(
    y = yc, X = Xc, 
    psi0 = as.numeric(coef(lm(yc ~ Xc - 1))), 
    prior_mean = prior_mean,
    prior_sd = prior_sd,
    burn = burn, 
    n_samples = n_samples,
    seed = 773
  )
}, error = function(e) {
  message("Warning: BETEL regression failed with error: ", e$message)
  return(NULL)
})

bhat_bt <- if (!is.null(fout_bt)) {
  mcmc_summary(fout_bt)
} else {
  message("Note: Using NA values for bhat_bt due to estimation failure.")
  matrix(NA, nrow = ncol(Xc), ncol = 4, 
         dimnames = list(colnames(Xc), c("est", "2.5 %", "97.5 %", "ESS")))
}


## =============================================================================
## plot
prepare_table <- function(df, var_names, method_name) {
  df_out <- as.data.frame(df) %>%
    setNames(c("est", "low", "upp", "ESS")) %>%
    mutate(Variable = var_names,
           Method = method_name) %>%
    select(Variable, Method, est, low, upp, ESS)
  
  return(df_out)
}

plot_data <- bind_rows(
  prepare_table(bhat_mle, cnames, "MLE"),
  prepare_table(bhat_gb, cnames, "gBayes"),
  prepare_table(bhat_bt, cnames, "BETEL"),
  prepare_table(bhat_cdp, cnames, "cDP"),
  prepare_table(bhat_cdp_nrm, cnames, "cDP-N")
)

plot_data$Variable <- factor(plot_data$Variable, 
                             levels = rev(unique(plot_data$Variable)))


legend_info <- plot_data %>%
  group_by(Method) %>%
  summarise(med_ess = median(ESS, na.rm = TRUE), .groups = 'drop') %>%
  mutate(new_label = paste0(Method, "(ESS: ", 
                            ifelse(is.na(med_ess), "NA", round(med_ess, 0)), ")"))

plot_data_final <- plot_data %>%
  left_join(legend_info, by = "Method") %>%
  mutate(Method_Label = factor(new_label, levels = legend_info$new_label))

ggplot(plot_data_final, aes(x = est, y = Variable, color = Method_Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  
  geom_errorbar(aes(xmin = low, xmax = upp), 
                width = 0.4, 
                linewidth = 0.8, 
                position = position_dodge(width = 0.7),
                orientation = "y") +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  
  theme_bw(base_size = 13) + 
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold", color = "black"),
    
    axis.title.x = element_text(size = 12, face = "bold", color = "black", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 1.2, color = "black"), 
    axis.ticks = element_line(linewidth = 1, color = "black")      
  ) +
  labs(
    x = "",
    y = ""
  ) +
  scale_color_brewer(palette = "Set1")
