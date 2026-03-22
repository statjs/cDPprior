
mle_regression <- function(y, X) {
  lmout <- lm(y ~ X - 1)
  cbind(est = coef(lmout), 
        confint(lmout))
}