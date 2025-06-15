merge_dists <- function(df1, df2, ratio) {
  
  atoms1 = df1$atoms
  atoms2 = df2$atoms
  weights1 = df1$weights
  weights2 = df2$weights
  
  adjusted_weights1 <- weights1 * ratio
  adjusted_weights2 <- weights2 * (1 - ratio)
  
  combined_atoms <- c(atoms1, atoms2)
  combined_weights <- c(adjusted_weights1, adjusted_weights2)
  
  combined_weights <- combined_weights / sum(combined_weights)
  
  return(data.frame(atoms = combined_atoms, weights = combined_weights))
}


rdp = function(total_mass, base_pr_measure, N=100, ... ) {
  
  betas = rbeta(N, 1, total_mass)
  weights = betas*cumprod(1-betas)
  weights[N] = 1 - sum(weights[1:(N-1)])
  atoms = base_pr_measure(N, ... )
  
  return(data.frame(weights=weights, atoms=atoms))
}


rbb = function(x) {

  n = length(x)
  atoms = x
  weights = rdirichlet(1, alpha=rep(1,n))[1,]
  df = data.frame(weights=weights, atoms=atoms)
  
  return(df)
}


rpostdp = function(x, total_mass, base_pr_measure, N=100, ... ) {
  
  df_prior = rdp(total_mass, base_pr_measure, N=100, ... )
  
  df_bb = rbb(x) 
  
  n = length(x)
  ratio = rbeta(1, n, total_mass)
  
  df = merge_dists(df_bb, df_prior, ratio)
  
  return(df)
}


mean_dist <- function(df) {
  sum(df$atoms * df$weights)
}

var_dist <- function(df) {
  m <- mean_dist(df)
  sum(df$weights * (df$atoms - m)^2)
}

sd_dist <- function(df) {
  v <- var_dist(df)
  sqrt(v)
}

skew_dist <- function(df) {
  m <- mean_dist(df)
  s = sd_dist(df)
  sum(df$weights * ((df$atoms - m)/s)^3)
}

kurt_dist <- function(df) {

  m <- mean_dist(df)
  s = sd_dist(df)
  sum(df$weights * ((df$atoms - m)/s)^4) 
}


quant_dist <- function(df, probs) {
  cdf <- data.frame(atoms = df$atoms, weights = cumsum(df$weights))
  sapply(probs, function(p) {
    atom <- cdf$atoms[which.max(cdf$weights >= p)]
    return(atom)
  })
}


rpost_bb_4moments = function(x, m) {
  
  mn = numeric(m)
  sig = numeric(m)
  skew = numeric(m)
  kurt = numeric(m) 
  
  for(i in 1:m) {
    dfn = rbb(x)
    mn[i] = mean_dist(dfn)
    sig[i] = sd_dist(dfn)
    skew[i] = skew_dist(dfn)
    kurt[i] = kurt_dist(dfn)
  }
  
  return(data.frame(mn = mn, sig=sig, skew=skew, kurt=kurt))
}


rpost_dp_4moments = function(x, m, total_mass=1, base_pr_measure=rcauchy, N=100, ... ) {
  
  mn = numeric(m)
  sig = numeric(m)
  skew = numeric(m)
  kurt = numeric(m) 
  
  for(i in 1:m) {
    dfn = rpostdp(x, total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, ...)
    mn[i] = mean_dist(dfn)
    sig[i] = sd_dist(dfn)
    skew[i] = skew_dist(dfn)
    kurt[i] = kurt_dist(dfn)
  }
  
  return(data.frame(mn = mn, sig=sig, skew=skew, kurt=kurt)) 
}



rprior_dp_4moments=function(m=2000, total_mass=1, base_pr_measure=rcauchy, N=100,...) {
  
  mn = numeric(m)
  sig = numeric(m)
  skew = numeric(m)
  kurt = numeric(m) 
  
  for(i in 1:m) {
    dfn = rdp(total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, ...)
    mn[i] = mean_dist(dfn)
    sig[i] = sd_dist(dfn)
    skew[i] = skew_dist(dfn)
    kurt[i] = kurt_dist(dfn)
  }
  
  prior_dp= data.frame(mn = mn, sig=sig, skew=skew, kurt=kurt)
  
  return(prior_dp)
}



prior_density = function(x, prior_par) {
  
  mn1=prior_par$mn1
  sig1=prior_par$sig1 
  a2=prior_par$a2 
  b2=prior_par$b2 
  mn3=prior_par$mn3 
  sig3=prior_par$sig3 
  a4=prior_par$a4 
  b4=prior_par$b4 
  
  if(is.vector(x)) {
    res = dcauchy(x[1], location=mn1, scale=sig1)*dgamma(x[2], shape=a2, rate=b2)*dcauchy(x[3], location=mn3, scale=sig3)*dgamma(x[4], shape=a4, rate=b4)
  } else {
    k = dim(x)[1]
    res = numeric(k)
    for(i in 1:k) {
      res[i] = dcauchy(x[i,1], location=mn1, scale=sig1)*dgamma(x[i,2], shape=a2, rate=b2)*dcauchy(x[i,3], location=mn3, scale=sig3)*dgamma(x[i,4], shape=a4, rate=b4)
    }
  } 
  return(res)
  
}



rpost_cdp_4moments1 = function(x, m, total_mass=1, base_pr_measure=rcauchy, N=100, prior_par=NULL, ... ) {
  
  post_dp = rpost_dp_4moments(x=x,m=m,total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, ...)
  
  prior_dp = rprior_dp_4moments(m=m, total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, ...)
  prior_dp_density = kde(prior_dp)  # 이 부분이 시간이 좀 걸린다. 
  density_denominator = predict(prior_dp_density, x= post_dp )
  density_denominator = pmax(density_denominator, 1e-300)
  
  if (is.null(prior_par)) {
    mn1 = mean(post_dp[,1])
    sig1 = 10*sd(post_dp[,1])
    
    m2= mean(post_dp[,2])
    s2= 10*sd(post_dp[,2])
    a2 = m2^2/s2^2
    b2 = m2/s2^2
    
    mn3 = mean(post_dp[,3])
    sig3 = 10*sd(post_dp[,3])
    
    m4= mean(post_dp[,4])
    s4= 10*sd(post_dp[,4])
    a4 = m4^2/s4^2
    b4 = m4/s4^2
    
    prior_par = list(mn1=mn1, sig1=sig1, a2=a2, b2=b2, mn3=mn3, sig3=sig3, a4=a4, b4=b4)
  }
  
  density_numerator = prior_density(post_dp, prior_par=prior_par)
  
  weights = density_numerator/density_denominator
  weights[is.na(weights) | is.infinite(weights)] <- 0
  weights = pmax(weights, 1e-300)
  indices = sample(length(weights), size=m, prob=weights)
  
  return(post_dp[indices, ])
  
}



silverman_bw = function(x) {
  res = 0.9*min(sd(x), IQR(x)/1.34)*length(x)^{1/5}
  return(res)
}



rpost_cdp_4moments2 = function(x, m, total_mass=1, base_pr_measure=rcauchy, N=100, prior_par=NULL, ... ) {

  post_dp = rpost_dp_4moments(x=x,m=m,total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, ...)
  
  prior_dp = rprior_dp_4moments(m=m, total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, ...)
  
  prior_dp_density_mn = kde(prior_dp$mn, h=silverman_bw(prior_dp$mn))  
  density_denominator_mn = predict(prior_dp_density_mn, x= post_dp$mn )
  
  prior_dp_density_sig = kde(prior_dp$sig, h=silverman_bw(prior_dp$sig))  
  density_denominator_sig = predict(prior_dp_density_sig, x= post_dp$sig )
  
  prior_dp_density_skew = kde(prior_dp$skew, h=silverman_bw(prior_dp$sig))  
  density_denominator_skew = predict(prior_dp_density_skew, x= post_dp$skew )
  
  prior_dp_density_kurt = kde(prior_dp$kurt, h=silverman_bw(prior_dp$kurt))  
  density_denominator_kurt = predict(prior_dp_density_kurt, x= post_dp$kurt )
  
  if (is.null(prior_par)) {
    mn1 = mean(post_dp[,1])
    sig1 = 10*sd(post_dp[,1])
    
    m2= mean(post_dp[,2])
    s2= 10*sd(post_dp[,2])
    a2 = m2^2/s2^2
    b2 = m2/s2^2
    
    mn3 = mean(post_dp[,3])
    sig3 = 10*sd(post_dp[,3])
    
    m4= mean(post_dp[,4])
    s4= 10*sd(post_dp[,4])
    a4 = m4^2/s4^2
    b4 = m4/s4^2
    
    prior_par = list(mn1=mn1, sig1=sig1, a2=a2, b2=b2, mn3=mn3, sig3=sig3, a4=a4, b4=b4)
  }
  
  
  density_numerator_mn = dcauchy(post_dp$mn, location=prior_par$mn1, scale=prior_par$sig1)
  density_numerator_sig = dgamma(post_dp$sig, shape=prior_par$a2, rate=prior_par$b2)
  density_numerator_skew = dcauchy(post_dp$skew, location=prior_par$mn3, scale=prior_par$sig3)
  density_numerator_kurt = dgamma(post_dp$kurt, shape=prior_par$a4, rate=prior_par$b4)
  
  weights_mn = density_numerator_mn/density_denominator_mn
  indices_mn = sample(length(weights_mn), size=m, prob=weights_mn)
  post_dp2_mn = post_dp$mn[indices_mn]
  
  weights_sig = density_numerator_sig/density_denominator_sig
  indices_sig = sample(length(weights_sig), size=m, prob=weights_sig)
  post_dp2_sig = post_dp$sig[indices_sig]
  
  weights_skew = density_numerator_skew/density_denominator_skew
  indices_skew = sample(length(weights_skew), size=m, prob=weights_skew)
  post_dp2_skew = post_dp$skew[indices_skew]
  
  weights_kurt = density_numerator_kurt/density_denominator_kurt
  indices_kurt = sample(length(weights_kurt), size=m, prob=weights_kurt)
  post_dp2_kurt = post_dp$kurt[indices_kurt]
  
  return(data.frame(mn=post_dp2_mn, sig=post_dp2_sig, skew=post_dp2_skew, kurt=post_dp2_kurt ))
  
}



rpost_cdp_4moments = function(x, m, sep=TRUE, total_mass=1, base_pr_measure=rcauchy, N=100, prior_par=NULL, ... ) {
  
  if(sep==FALSE) {
    return(rpost_cdp_4moments1(x=x, m=m, total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, prior_par=prior_par, ... ))
  } else {
    return(rpost_cdp_4moments2(x=x, m=m, total_mass=total_mass, base_pr_measure=base_pr_measure, N=N, prior_par=prior_par, ... ))
  }
}


