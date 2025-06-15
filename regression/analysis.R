source("mcmc.R")

# 패키지 로드
library(GGally)
library(ggplot2)
library(tidyverse)
library(ggcorrplot)
library(mvtnorm) # rmvnorm(), rmvt() 
library(gtools) # rdirichlet() 
library(coda)
library(ggmcmc)
library(ks) #kde() 

### 데이터
load("nhanes.RData")


### 산점도 행렬 시각화
plot_data <- nhanes_clean %>%
  select(
    LBXGLU,       # 공복 혈당
    BMXWT,        # 체중
    BMXWAIST,     # 허리둘레
    BPXOSY1,      # 수축기 혈압
    BPXODI1,      # 이완기 혈압
    LBXIN,        # 인슐린
    RIDAGEYR,     # 나이
    PHAFSTHR      # 공복 시간
  )

# 산점도 행렬 시각화
ggpairs(plot_data)


# 상관행렬 히트맵
cor_mat <- cor(nhanes_clean[, -1])  # SEQN 제외
ggcorrplot(cor_mat, lab = TRUE)


### 회귀모형 적합
y = nhanes_clean %>% select(LBXGLU) %>% unlist()
Y = matrix(y, ncol=1) # 반응변수 행렬 n x 1
X = nhanes_clean %>% select(-c(LBXGLU, SEQN, BMXHT, RIAGENDR)) %>% as.matrix() # 설명변수 행렬 n x p 
n = dim(X)[1]
p = dim(X)[2]


library(mvcauchy)
set.seed(13)

total_mass=1
base_pr_measure=rmvcauchy
N=100
m=1000 
post_cdp = rpost_cdp_regcoef(Y=Y, X=X, m=1000, 
                             total_mass=1, base_pr_measure=base_pr_measure, N=N, 
                             mu=rep(0,p+1), sigma=10*diag(p+1), 
                             prior_par=list(mean=rep(0,p), sig=rep(10,p))) 
post_cdp %>% mcmc %>% summary


### 빈도론적 방법
# 샘플 공분산 계산
S_xx <- cov(X)
S_xy <- cov(X, y)

# 회귀계수 추정
beta_hat <- solve(S_xx, S_xy)
beta_hat

# Bootstrap 신뢰구간
set.seed(123)
B <- 1000
n <- nrow(X)
p <- ncol(X)

beta_boot <- matrix(NA, nrow = B, ncol = p)

for (b in 1:B) {
  idx <- sample(1:n, size = n, replace = TRUE)
  X_b <- X[idx, ]
  y_b <- y[idx]
  
  S_xx_b <- cov(X_b)
  S_xy_b <- cov(X_b, y_b)
  
  beta_boot[b, ] <- solve(S_xx_b, S_xy_b)
}

# 신뢰구간 계산 (percentile method)
conf_int <- apply(beta_boot, 2, quantile, probs = c(0.025, 0.975))
rownames(conf_int) <- c("2.5%", "97.5%")
colnames(conf_int) <- colnames(X)
conf_int


# 테이블 데이터 프레임 생성
beta_df <- data.frame(
  Variable = colnames(X),
  Estimate = round(beta_hat, 4),
  `2.5%` = round(conf_int[1, ], 4),
  `97.5%` = round(conf_int[2, ], 4)
)
beta_df



