library(tidyverse) # ggplot2, dplyr 패키지들의 모임 
library(GGally) ## ggpairs
library(Hmisc) # hist를 데이터프레임에 적용할 수 있게 해준다. 이산형자료 요약에 좋은 describe가 있다. 
library(psych) # 연속형 자료 요약에 좋은 describe가 있다. 
library(summarytools) #desc(), stby()
library(gtools) # rdirichlet() 
library(readxl)
library(coda)
library(ggmcmc)
library(ks) #kde() 

set.seed(1234567)

rm(list=ls())

somites <- c(79, 89, 93, 95, 99, 100, 102, 104, 105, 106, 107, 113, 114, 115, 118, 120, 122, 123, 124, 125,
             126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144,
             145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 159, 162, 164)

frequencies <- c(1, 2, 1, 1, 2, 1, 2, 3, 1, 1, 2, 1, 3, 3, 2, 5, 1, 2, 2, 5,
                 2, 1, 5, 4, 1, 4, 7, 4, 3, 5, 6, 9, 10, 6, 5, 17, 27, 21, 24,
                 37, 32, 31, 29, 26, 51, 15, 15, 12, 9, 9, 10, 4, 2, 2, 1)

n_total <- sum(frequencies)  

somite_data <- rep(somites, frequencies)

df <- data.frame(somite = somite_data)

ggplot(df, aes(x = somite)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(x = "Number of Somites", y = "Frequency") +
  theme_minimal()


### bootstrap
library(boot)
library(e1071)  # skewness 

set.seed(1)

# 1. 평균에 대한 부트스트랩
boot_mean <- boot(data = somite_data, statistic = function(x, i) mean(x[i]), R = 1000)
boot.ci(boot_mean, type = "perc")

# 2. 표준편차
boot_sd <- boot(data = somite_data, statistic = function(x, i) sd(x[i]), R = 1000)
boot.ci(boot_sd, type = "perc")

# 3. 왜도
boot_skew <- boot(data = somite_data, statistic = function(x, i) skewness(x[i]), R = 1000)
boot.ci(boot_skew, type = "perc")

# 4. 첨도
boot_kurt <- boot(data = somite_data, statistic = function(x, i) kurtosis(x[i]), R = 1000)
boot.ci(boot_kurt, type = "perc")

### 사후분포
source("mcmc.R")

compute_hyperparameters <- function(y, B = 1000) {
  n <- length(y)
  
  # 1. m1 and tau1
  m1 <- mean(y)
  
  # Bootstrap to get sd of sample median
  means <- replicate(B, mean(sample(y, replace = TRUE)))
  tau1 <- 10 * sd(means) * 2
  
  # 2. a2 and b2 
  gamma_hat <- sd(y)
  
  # Bootstrap to get variance of IQR/2
  gamma_hats <- replicate(B, sd(sample(y, replace = TRUE)))
  var_gamma <- var(gamma_hats)
  
  a2 <- (gamma_hat^2) / (100 * var_gamma * 4)
  b2 <- gamma_hat / (100 * var_gamma * 4)
  
  # 3. m3 and tau3 
  y_std <- (y - m1) / gamma_hat
  m3 <- mean(y_std^3)
  
  # Bootstrap sd of skewness
  m3_boot <- replicate(B, {
    y_b <- sample(y, replace = TRUE)
    y_std_b <- (y_b - mean(y_b)) / (sd(y_b))
    mean(y_std_b^3)
  })
  tau3 <- 10 * sd(m3_boot) * 2
  
  # 4. a4 and b4 
  m4 <- mean(y_std^4)
  
  m4_boot <- replicate(B, {
    y_b <- sample(y, replace = TRUE)
    y_std_b <- (y_b - mean(y_b)) / (sd(y_b))
    mean(y_std_b^4)
  })
  var_m4 <- var(m4_boot)
  
  a4 <- (m4^2) / (100 * var_m4 * 4)
  b4 <- m4 / (100 * var_m4 * 4)
  
  list(
    mn1 = m1, sig1 = tau1,
    a2 = a2, b2 = b2,
    mn3 = m3, sig3 = tau3,
    a4 = a4, b4 = b4
  )
}


prior_par = compute_hyperparameters(y=somite_data)

library(extraDistr)

set.seed(1377)
post_cdp = rpost_cdp_4moments(x=somite_data, m=1000, sep=TRUE, total_mass=1,
                              base_pr_measure=rlst, df=5, mu=median(somite_data), sigma=IQR(somite_data)/1.34, N=100,
                              prior_par = prior_par)
post_cdp %>% mcmc %>% summary

