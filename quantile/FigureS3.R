library(ggplot2)
source("dp.R")

rt_iqr <- function(n, df = 3, target_iqr = 9) {
  s <- target_iqr / (qt(0.75, df) - qt(0.25, df))
  s * rt(n, df = df)
}

n_samples <- 10000
burn <- 5000

B <- 100
n_grid <- c(100, 500)
target_iqr <- 9
df <- 3

res_list <- vector("list", length(n_grid) * B)
k <- 1

set.seed(1)

for (n in n_grid) {
  for (b in 1:B) {
    
    # --- data
    y <- rt_iqr(n = n, df = df, target_iqr = target_iqr)
    
    # --- cDP (IQR constraint)
    fout_cDP <- cDPprior::cdp_quantile(x = y, p = 0.5, n_samples = n_samples, base_dist = 1, 
                                       iqr = target_iqr)
    med_cdp <- mean(fout_cDP$samples[(burn+1):n_samples])   
    
    # --- DP baseline
    fout_DP <- dp_post_median_iqr(x = y, n_samples = n_samples, alpha = 1, K = 100)
    fout_DP <- fout_DP[((burn+1):n_samples),,drop=FALSE]
    med_dp  <- mean(fout_DP$q50, na.rm = TRUE)
    iqr_dp <- mean(fout_DP$iqr,  na.rm = TRUE)
    
    # --- store (posterior means only)
    res_list[[k]] <- data.frame(
      n = n, rep = b,
      med_cdp = med_cdp,
      iqr_cdp = target_iqr,
      med_dp = med_dp,
      iqr_dp = iqr_dp
    )
    k <- k + 1
  }
  message("done n = ", n)
}

res <- do.call(rbind, res_list)
aggregate(cbind(med_cdp, med_dp, iqr_dp) ~ n, data = res,
          FUN = function(z) c(mean = mean(z, na.rm=TRUE), sd = sd(z, na.rm=TRUE)))
res$n <- factor(res$n, levels = c(100, 500))

## --- mean comparison 
res100 <- subset(res, n == 100)
df <- rbind(
  data.frame(method = "DP",  med = res100$med_dp,  iqr = res100$iqr_dp),
  data.frame(method = "cDP", med = res100$med_cdp, iqr = 9)
)
df$method <- factor(df$method, levels = c("cDP","DP"))

paper_theme <- theme_classic(base_size = 14) +
  theme(
    text = element_text(color="black"),
    axis.title = element_text(face="bold"),
    axis.text  = element_text(color="black"),
    plot.title = element_text(face="bold", hjust=0.5),
    axis.line  = element_line(linewidth=0.9),
    axis.ticks = element_line(linewidth=0.9),
    legend.position = "none"
  )

jrss_cols <- c(
  "cDP" = "#6E6E6E",  # medium grey (더 밝게)
  "DP"  = "#C7C7C7"   # light grey
)

# --- Mean
ggplot(df, aes(x = method, y = med, fill = method)) +
  geom_boxplot(width = 0.6, linewidth = 0.8, outlier.size = 0.8) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.9) +
  labs(x = NULL, y = expression("Posterior mean of " ~ Q[0.5](F))) +
  scale_fill_manual(values = jrss_cols) +
  paper_theme


# --- IQR (DP only, same color)
ggplot(subset(df, method == "DP"),
       aes(x = method, y = iqr, fill = method)) +
  geom_boxplot(width = 0.6, linewidth = 0.8, outlier.size = 0.8) +
  geom_hline(yintercept = 9, linetype = 2, linewidth = 0.9) +
  labs(x = NULL, y = "Posterior mean of IQR") +
  scale_fill_manual(values = jrss_cols) +
  paper_theme

save.image("median_iqr_dp.RData")


