# Find quantile value for normal mixture at quantile level(q).
# q: quantile level
# mu: means of normal distributions for normal mixture distribution
# sd: standard deviations of normal distributions for normal mixture distribution
# pi: weights for normal distributions for normal mixture distribution

FmixNorm=function(q,mu,sd,pi){
  mu=as.vector(mu)
  sd=as.vector(sd)
  pi=as.vector(pi)
  qres=NULL
  for(i in 1:length(mu))
  {
    qres = c(qres, pi[i]*pnorm(q, mean = mu[i], sd = sd[i]))
  }
  return(sum(qres))
}

QmixNorm = function(q,mu,sd,pi){
  minmax <- range(qnorm(q,mu[1],sd[1]),qnorm(q,mu[2],sd[2]), qnorm(q,mu[3],sd[3])) 
  uniroot(function(x) FmixNorm(x,mu,sd,pi)-q,
          interval = minmax,
          tol = 10^{-16})$root  
}
 