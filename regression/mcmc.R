reg_coef = function(Y, X, weights, as_vector=FALSE) {
  
  X = as.matrix(X) 
  n = dim(X)[1]
  p = dim(X)[2]
  A1 = matrix(0, nrow=p, ncol=p) 
  A2 = matrix(0, nrow=p, 1)
  for(i in 1:n) {
    xtemp = as.matrix(X[i,])
    A1 = A1 + weights[i]*xtemp %*% t(xtemp) 
    A2 = A2 + weights[i]*xtemp*Y[i,1]
  }
  beta = solve(A1) %*% A2
  
  if(as_vector == TRUE) {
    return(as.vector(beta))
  } else {
    return(beta)
  }
}

regcoef_dist <- function(distn, as_vector=FALSE) {
  return(reg_coef(Y=distn$Y, X=distn$X, weights=distn$weights, as_vector=as_vector))
}


rpost_bb_regcoef = function(Y, X, m=1000) {
  
  n = dim(X)[1]
  p = dim(X)[2]
  post_beta = matrix(0, nrow=m, ncol=p)
  for(i in 1:m) {
    weights = rdirichlet(1, alpha=rep(1,n))[1,]
    post_beta[i,] = reg_coef(Y=Y, X=X, weights=weights)
  }
  return(data.frame(post_beta))
}



rbb_reg = function(Y, X) {
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  weights = rdirichlet(1, alpha=rep(1,n))[1,]
  distn = list(weights=weights, Y=Y, X=X)
  
  return(distn)
}


rdp_reg = function(total_mass, base_pr_measure, N=100, ... ) {

  betas = rbeta(N, 1, total_mass)
  weights = betas*cumprod(1-betas)
  weights[N] = 1 - sum(weights[1:(N-1)])
  atoms = base_pr_measure(N, ... )
  p = dim(atoms)[2]-1
  Y = matrix(atoms[,1], ncol=1)
  X = atoms[, 2:(p+1)]
  
  return(list(weights=weights, Y=Y, X=X))
}


merge_distns_reg <- function(distn1, distn2, ratio) {

  combined_Y <- rbind(as.matrix(distn1$Y), as.matrix(distn2$Y))
  combined_X <- rbind(as.matrix(distn1$X), as.matrix(distn2$X))
  combined_weights <- c(distn1$weights*ratio,distn2$weights*(1-ratio))
  
  # 결과 정규화 (합이 1이 되도록)
  combined_weights <- combined_weights / sum(combined_weights)
  
  return(list(weights = combined_weights, Y=combined_Y, X=combined_X))
}


rpostdp_reg = function(Y, X, total_mass, base_pr_measure, N=100, ... ) {
  
  distn_prior = rdp_reg(total_mass = total_mass, base_pr_measure=base_pr_measure, N=N, ... )
  
  distn_bb = rbb_reg(Y=Y, X=X) 
  
  n = dim(X)[1]
  ratio = rbeta(1, n, total_mass)
  
  distn = merge_distns_reg(distn_bb, distn_prior, ratio)
  
  return(distn)
}


rpost_dp_regcoef = function(Y, X, m, total_mass=1, base_pr_measure=rmvnorm, N=100, ... ) {
  
  n = dim(X)[1]
  p = dim(X)[2]
  post_beta = matrix(0, nrow=m, ncol=p)
  for(i in 1:m) {
    distn = rpostdp_reg(Y=Y,X=X, total_mass = total_mass, base_pr_measure = base_pr_measure, N=N, ...)
    post_beta[i,] = regcoef_dist(distn, as_vector=TRUE)
  }
  return(data.frame(post_beta))
  
}


rprior_dp_regcoef = function(m=1000, total_mass, base_pr_measure, N=100, ...) {
  
  distn = rdp_reg(total_mass = total_mass, base_pr_measure=base_pr_measure, N=N, ...)
  beta = regcoef_dist(distn, as_vector=TRUE)
  p = length(beta)
  
  prior_beta = matrix(0,nrow=m, ncol=p)
  
  for(i in 1:m) {
    distn = rdp_reg(total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, ...)
    prior_beta[i,] = regcoef_dist(distn, as_vector=TRUE)
  }
  
  return(data.frame(prior_beta))
}


silverman_bw = function(x) {
  res = 0.9*min(sd(x), IQR(x)/1.34)*length(x)^{1/5}
  return(res)
}


rpost_cdp_regcoef = function(Y,X, m, total_mass=1, base_pr_measure=rmvnorm, N=100, prior_par=NULL, ... ) {

  n = dim(X)[1]
  p = dim(X)[2]
  
  post_dp = rpost_dp_regcoef(Y=Y,X=X,m=m,total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, ...)
  
  prior_dp = rprior_dp_regcoef(m=1000, total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, ...)
  
  post_beta = matrix(0,nrow=m, ncol=p)
  for(j in 1:p) {
    prior_dp_density = kde(prior_dp[,j], h=silverman_bw(prior_dp[,j]))  
    density_denominator = predict(prior_dp_density, x= post_dp[,j] )
    if(is.null(prior_par)) {
      weights = 1/density_denominator
    } else {
      density_numerator = dcauchy(post_dp[,j], location=prior_par$mean[j], scale=prior_par$sig[j])
      weights = density_numerator/density_denominator
    }
    indices = sample(length(weights), size=m, prob=weights)
    post_beta[,j] = post_dp[indices,j]
  }
  
  
  return(data.frame(post_beta) )
  
}


