library(tidyverse)

hurricane <- readr::read_csv("HourlyMaxWinds18512006.csv")
hurricane

hurricane %>% select(Yr) %>% unique() %>% range()

# histogram for data
hurricane %>% select(Yr, Wmax) %>% subset(Yr == 2006) %>% 
  ggplot(aes(x = Wmax)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "skyblue", color = "black") +
  labs(x = "Maximum Wind Speed (Wmax)", y = "Density") +
  theme_minimal() + theme_bw()
ggsave("./figs/wmax_hist_2006.png", width = 8, height = 5, dpi = 300)



# sample quantile
q_vals <- hurricane %>%
  filter(Yr == 2006) %>%
  summarise(q95 = quantile(Wmax, 0.95),
            q99 = quantile(Wmax, 0.99)) %>%
  pivot_longer(everything(), names_to = "quantile", values_to = "value") %>%
  mutate(label = ifelse(quantile == "q95", "95% Quantile", "99% Quantile"),
         color = ifelse(quantile == "q95", "orange", "red"))

### Slice sampler
source("mcmc.R")

x <- hurricane %>%
  filter(Yr == 2006) %>% select(Wmax) %>% 
  unlist()

# 95%, 99%
n_samples <- 10000
samples <- list()
p_seq <- c(0.95, 0.99)
for (i in 1:length(p_seq)) {
  samples[[i]] <- slice_sampler(
    log_posterior = log_posterior_cauchy, # log_posterior_nrm, log_posterior_cauchy
    theta0 = quantile(x, p_seq[i]),
    w = 0.5,
    n_samples = n_samples,
    max_step_out = 10,
    x = x,
    mu0 = median(x), sigma0 = IQR(x) / 2,
    m = quantile(x, p_seq[i]), tau = 10*(pi*(IQR(x)/2)/sqrt(length(x)))*(1 + tan(pi*(p_seq[i]-0.5))^2)*1.349,
    p = p_seq[i]
  )
}
n_burnin <- ceiling(n_samples/2)
post_theta <- data.frame(p_0_95 = samples[[1]][(n_burnin+1):n_samples], p_0_99 = samples[[2]][(n_burnin+1):n_samples])


# mcmc plot
library(ggmcmc)
library(coda)

post_theta %>% mcmc %>% ggs %>% ggs_traceplot() +
  ylab(expression(theta))
# ggsave("./figs/posterior_traceplot.png", width = 8, height = 5, dpi = 300)


post_theta %>% mcmc %>% ggs %>% ggs_density() +
  ylab(expression(theta))
# ggsave("./figs/posterior_density.png", width = 8, height = 5, dpi = 300)


post_theta %>% mcmc %>% ggs %>% ggs_pairs() +
  ylab(expression(theta))


# posterior summary
post_theta %>% mcmc %>% summary()

q_post <- post_theta %>%
  summarise(
    P_0_95_mean = mean(p_0_95),
    P_0_95_lower = quantile(p_0_95, 0.025),
    P_0_95_upper = quantile(p_0_95, 0.975),
    P_0_99_mean = mean(p_0_99),
    P_0_99_lower = quantile(p_0_99, 0.025),
    P_0_99_upper = quantile(p_0_99, 0.975)
  ) %>%
  pivot_longer(everything(), names_to = c("quantile", "stat"), names_sep = "_", names_prefix = "P_0_") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(color = ifelse(quantile == "95", "orange", "red"),
         label = ifelse(quantile == "95", "95% Quantile", "99% Quantile"))



hurricane %>%
  filter(Yr == 2006) %>%
  ggplot(aes(x = Wmax)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, fill = "skyblue", color = "black") +
  
  # Credible interval (shaded area)
  geom_rect(data = q_post,
            aes(xmin = lower, xmax = upper, ymin = 0, ymax = Inf, fill = label),
            alpha = 0.2, inherit.aes = FALSE) +
  
  # Posterior mean
  geom_vline(data = q_post,
             aes(xintercept = mean, color = label),
             linetype = "dashed", linewidth = 1) +
  
  # Posterior mean label
  geom_text(data = q_post,
            aes(x = mean, y = 0, label = round(mean, 1), color = label),
            vjust = 1.5, hjust = 0.5, size = 4, show.legend = FALSE) +
  
  # Sample quantile (black dashed line)
  geom_vline(data = q_vals,
             aes(xintercept = value),
             color = "black", linetype = "dotted", linewidth = 1) +
  
  scale_color_manual(name = "Posterior Mean",
                     values = c("95% Quantile" = "orange", "99% Quantile" = "red")) +
  scale_fill_manual(name = "95% Credible Interval",
                    values = c("95% Quantile" = "orange", "99% Quantile" = "red")) +
  
  labs(x = "Maximum Wind Speed (Wmax)",
       y = "Density") +
  theme_bw()
# ggsave("./figs/posterior_cauchy.png", width = 8, height = 5, dpi = 300)

q_post