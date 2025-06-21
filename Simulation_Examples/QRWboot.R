library(quantreg)
source('cv.pred.qsvcm.R')

#' Calculate the confidence intervals for quantile spatially varying coefficients using the BPST.
#'
#' \code{QRWboot} constructs the confidence intervals of quantile spatially varying coefficients.
#'
#' @param y The response of dimension \code{n} by one, where \code{n} is the number of observations.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param S The cooridinates of dimension \code{n} by two. Each row is the coordinates of an observation.
#' \cr
#' @param V The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 2.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param BQ2_pop Population bivariate spline basis times Q2 matrix in the paper, that is B.pop%*%Q2. 
#' \cr
#' @param lambda The vector of the candidates of penalty parameter.
#' \cr
#' @param tau The quantile level -- default is 0.5 (median regression).
#' \cr
#' @param eta The tuning (penalty) parameter for each iteration in the ADMM ($\eta^{(j)}$ in the manuscript) -- default is 1.0.
#' \cr
#' @param eps.abs The absolute tolerances in the ADMM -- default is 10^(-4). 
#' \cr
#' @param eps.rel The relative tolerances in the ADMM -- default is 10^(-2).
#' \cr
#' @param step Type of varying method for the ADMM. There are non-varying and varying methods (1: non-varying vs 2: varying) -- default is 1.
#' \cr
#' @param mu.eta Varying tuning parameter related to the ADMM ($\zeta^{\ast}$ in the manuscript) -- default is 10. 
#' \cr
#' @param incr Varying tuning parameter related to the ADMM ($\delta^{\incr}$ in the manuscript) -- default is 2. 
#' \cr
#' @param B The bivariate spline basis -- default is NULL.
#' \cr
#' @param Q2 The matrix Q2 for the reparametrization of the spline coefficients in the paper -- default is NULL.
#' \cr
#' @param P The block diagonal penalty matrix P in the paper -- default is NULL. 
#' \cr
#' @param nB The number of bootstrap iterations -- default is 100.
#' \cr
#' 
#' @return The function returns a list with the following items:
#' \item{Lower.CB}{The lower bound of confidence intervals of quantile spatially varying coefficients on the domain.}
#' \item{Upper.CB}{The upper bound of confidence intervals of quantile spatially varying coefficients on the domain.}



QRWboot = function(y, X, S, V, Tr, d = 2, r = 1, BQ2_pop, lambda = 10^(seq(log10(10^(-6)), log10(10^(1)), length.out = 10)), 
                   tau = 0.5, eta = 1.0, 
                   eps.abs = 10^(-4), eps.rel = 10^(-2), step = 1, mu.eta = 10, incr = 2,
                   B = NULL, Q2 = NULL, P= NULL, nB = 100)  
{
  
  
  
  lambda.cv = lambda
  all.cv.result = list()
  all.cv.measure = c()
  for(lam_iter in 1:nlambda){
    cvlam = lambda.cv[lam_iter] 
    cvresult = cv.pred.qsvcm(y, X, S, V, Tr, d = d, r = r, lambda = cvlam, tau = tau, nfold = 5)
    all.cv.result[[lam_iter]] = cvresult
    cv.measure = c(mean(cvresult[[1]]), mean(cvresult[[2]]))
    all.cv.measure =  rbind(all.cv.measure, cv.measure)
  }
  
  
  # set up for smoothing parameters in the penalty term (update)
  initial.cv.result = cbind(round(lambda.cv, 2), all.cv.measure)
  ind.cv = which.min(initial.cv.result[,2])
  initial.lambda = lambda.cv[ind.cv]
  
  new.nlambda = 10 
  lambda.scale = 100
  new.lambda = 10^(seq(log10(initial.lambda/lambda.scale),log10(initial.lambda*lambda.scale),length.out=new.nlambda))
  
  lambda.cv = new.lambda
  all.cv.result2 = list()
  all.cv.measure2 = c()
  for(lam_iter in 1:new.nlambda){
    cvlam = lambda.cv[lam_iter] 
    cvresult = cv.pred.qsvcm(y, X, S, V, Tr, d = d, r = r, lambda = cvlam, tau = tau, nfold = 5)
    all.cv.result2[[lam_iter]] = cvresult
    cv.measure = c(mean(cvresult[[1]]), mean(cvresult[[2]]))
    all.cv.measure2 =  rbind(all.cv.measure2, cv.measure)
  }
  
  final.cv.result = cbind(round(lambda.cv, 2), all.cv.measure2)
  ind.cv = which.min(final.cv.result[,2])
  final.lambda = lambda.cv[ind.cv]
  
  #
  lambda = matrix(final.lambda, 1)
  Reference.result = fit.qsvcm(y, X, S, V, Tr, d, r, lambda, tau, eta, eps.abs, eps.rel, step, mu.eta, incr, B, Q2, P)
  gamma = Reference.result$gamma  
  
  B = Reference.result$B
  Q2 = Reference.result$Q2
  P = Reference.result$P
  
  BQ2 = as.matrix(B%*%Q2)
  np = ncol(X)
  J = dim(BQ2)[2]
  
  W = as.matrix(kr(X,BQ2,byrow=TRUE))
  
  star_beta_all_all=list(B)
  star_pop_beta_all_all=list(B)
  
  for(biter in 1:nB){

    yhat = W%*%gamma
    
    beta_all = c()
    pop_beta_all = c()
    for(p in 1:np){
      thetamat = matrix(gamma, ncol=np) 
      
      beta = BQ2%*%thetamat[,p]
      beta_pop = BQ2_pop%*%thetamat[,p]
      
      beta_all = cbind(beta_all, beta)
      pop_beta_all = cbind(pop_beta_all, beta_pop)
    }
    tilte_beta = beta_all
    pop_tilte_beta = pop_beta_all

    hat_eps = y-yhat
    hat_eps = hat_eps + (1/akj(hat_eps,0)$dens)*diag(W%*%solve(t(W)%*%W+diag(lambda[1],dim(W)[2]))%*%t(W))*(tau-(hat_eps<0))
  
    set.seed(biter)
    r = c(2*(1-tau), -2*tau)
    pr = c(1-tau, tau)
    ri = sample(x = r, size = length(y), replace = TRUE, prob = pr)
    star_eps = ri * abs(hat_eps)

    star_y = yhat + star_eps

    boot.result = fit.qsvcm(star_y, X, S, V, Tr, d, r, lambda, tau, eta, eps.abs, eps.rel, step, mu.eta, incr, B, Q2, P)

    star_beta_all = c()
    star_pop_beta_all = c()
    for(p in 1:np){
      star_thetamat = matrix(boot.result$gamma, ncol=np) 
      star_beta = BQ2%*%star_thetamat[,p]
      star_pop_beta = BQ2_pop%*%star_thetamat[,p]
      
      star_beta_all = cbind(star_beta_all, star_beta)
      star_pop_beta_all = cbind(star_pop_beta_all, star_pop_beta)
    }
    star_tilte_beta = star_beta_all
    star_pop_tilte_beta = star_pop_beta_all
    
    star_beta_all_all[[biter]] = star_tilte_beta 
    star_pop_beta_all_all[[biter]] = star_pop_tilte_beta 
  }
  
  result_beta_all =c() 
  result_pop_beta_all =c() 
  
  for(biter in 1:nB){
    result_beta_all = cbind(result_beta_all, as.vector(star_beta_all_all[[biter]]))
    result_pop_beta_all = cbind(result_pop_beta_all, as.vector(star_pop_beta_all_all[[biter]]))
  }
  

  

  dist.beta = sqrt(length(y))*(result_beta_all - as.vector(tilte_beta))
  dist.pop.beta = sqrt(length(y))*(result_pop_beta_all - as.vector(pop_tilte_beta))

  
  lbound = apply(dist.beta, 1, quantile, probs = c(0.025))
  ubound = apply(dist.beta, 1, quantile, probs = c(0.975))
  
  pop.lbound = apply(dist.pop.beta, 1, quantile, probs = c(0.025))
  pop.ubound = apply(dist.pop.beta, 1, quantile, probs = c(0.975))
  
  
  BLB = as.vector(tilte_beta)-1/sqrt(length(y))*ubound
  BUB = as.vector(tilte_beta)-1/sqrt(length(y))*lbound
  SLower.CB = matrix(BLB, ncol=np)
  SUpper.CB = matrix(BUB, ncol=np)

  pop.BLB = as.vector(pop_tilte_beta)-1/sqrt(length(y))*pop.ubound
  pop.BUB = as.vector(pop_tilte_beta)-1/sqrt(length(y))*pop.lbound
  
  
  Lower.CB = matrix(pop.BLB, ncol=np)
  Upper.CB = matrix(pop.BUB, ncol=np)


  list(Lower.CB = Lower.CB, Upper.CB = Upper.CB) 
}
