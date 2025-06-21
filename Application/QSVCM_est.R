#' Calculate the confidence intervals for quantile spatially varying coefficients using the BPST.
#'
#' \code{QSVCM_est} constructs the confidence intervals of quantile spatially varying coefficients.
#'
#' @param y The response of dimension \code{n} by one, where \code{n} is the number of observations.
#' \cr
#' @param X The design matrix of dimension \code{n} by \code{p}, with an intercept. Each row is an observation vector.
#' \cr
#' @param B The bivariate spline basis.
#' \cr
#' @param Q2 The matrix Q2 for the reparametrization of the spline coefficients in the paper.
#' \cr
#' @param P The block diagonal penalty matrix P in the paper. 
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
#' @param method Type of algorithm for estimation -- default is "ADMM". The other choice is "ORADMM" called over-relaxed ADMM.
#' \cr
#' @param step Type of varying method for the ADMM. There are non-varying and varying methods (1: non-varying vs 2: varying) -- default is 1.
#' \cr
#' @param mu.eta Varying tuning parameter related to the ADMM ($zeta^{\ast}$ in the manuscript) -- default is 10. 
#' \cr
#' @param incr Varying tuning parameter related to the ADMM ($delta^{\incr}$ in the manuscript) -- default is 2. 
#' \cr
#' @param alpha A tuning parameter related to the ORADMM -- default is 1.5. 
#' \cr
#' @param b.initial Initial value for spline coefficient (gamma) for ADMM algorithm. 
#' \cr
#' 
#' @return The function returns a list with the following items:
#' \item{gamma}{Spline coefficients before multiplying the matrix Q2.}
#' \item{theta}{Spline coefficient after multiplying the matrix Q2. That is, Q2*gamma.}
#' \item{beta}{The estimated quantile spatially varying coefficient functions.}
#' \item{gacv}{The resulting value of criterion for smoothness parameter (GACV).}
#' \item{lambdac}{Chosen lambda based on criterion for smoothness parameter.}
#' \item{time}{Computational time.}
#' \item{iter}{The number of iterations for the ADMM.}
 

QSVCM_est=function(y,X,B,Q2,P,lambda, tau = 0.5, eta = 1.0, 
                   eps.abs = 10^(-4), eps.rel = 10^(-2),
                   method = "ADMM", step = 1, 
                   mu.eta = 10, incr = 2, alpha = 1.5,    
                   b.initial = FALSE
                   ){
  start_time = Sys.time()
  
  if(tau < 0 | tau > 1 ){
    warning("Quantile level should be between 0 and 1. The default quantile level of 0.5 is used.")
    tau = 0.5
  }

  np = ncol(X)
  J = dim(Q2)[2]

  lambda=as.matrix(lambda)
  nl=nrow(lambda)
  if(ncol(lambda)==1){
    lambda=matrix(rep(lambda,times=np),nl,np)
  }
  
 
  BQ2 = B%*%Q2
  W = as.matrix(kr(X,BQ2,byrow=TRUE))
  

  alpha_all = matrix(rep(0,(np*J)*nl),ncol=nl)
  gacv_all = rep(0,nl)
  df_all = rep(0,nl)
  eta_all = rep(0,nl)
  

  for(il in 1:nl){
    var.eta = eta
     
    a=1/(2*eta)

     solve.mat = solve((2/var.eta)*kronecker(diag(lambda[il,]),P)+crossprod(W)) 
     
    old.u=new.u=(1/eta)*y
    if(b.initial==FALSE){
 
      old.b=new.b = solve.mat%*%t(W)%*%y         

    } else {
      old.b=new.b=b.initial
    }
    v=-old.u+y-W%*%old.b-(2*tau-1)/(2*eta)
    old.r=new.r=stoper(v,a)

    iter.admm = 0
    all.var.eta=c()
    
    repeat{
      iter.admm = iter.admm + 1

      old.u=new.u
      old.b=new.b
      old.r=new.r

      v=-old.u+y-W%*%old.b-(2*tau-1)/(2*var.eta)
      new.r=stoper(v,a)
      
      if(method == "ADMM"){
        new.b=solve.mat%*%t(W)%*%(y-new.r-old.u) 
        new.u=old.u+(new.r+W%*%new.b-y)
        
      } else if(method == "ORADMM"){
        new.b=solve.mat%*%t(W)%*%(y-(alpha*new.r-(1-alpha)*(W%*%old.b-y))-old.u)  
        new.u=old.u+(alpha*new.r-(1-alpha)*(W%*%old.b-y)+W%*%new.b-y)
      }

      
      r.prime=y-W%*%new.b-new.r
      r.dual=var.eta*W%*%(new.b-old.b) #
      eps.prime=sqrt(length(r.prime))*eps.abs+eps.rel*max(norm(new.r, type="2"),norm(W%*%new.b, type="2"),norm(y, type="2"))
      eps.dual=sqrt(length(r.dual))*eps.abs+eps.rel*norm(new.u,type="2")
      
      
      if(step == 1){
        var.eta = var.eta
        new.u = new.u
        a = 1/ (2 * var.eta)
      } else if(step == 2){
        if(norm(r.prime, type="2") > mu.eta * norm(r.dual, type="2")){
          var.eta = incr*var.eta
          new.u = new.u/incr
          a = 1/ (2 * var.eta)
        } else if (mu.eta * norm(r.prime,type="2") < norm(r.dual,type="2")){
          var.eta = var.eta / incr
          new.u = incr*new.u
          a = 1/ (2 * var.eta)
        } else {
          var.eta = var.eta
          new.u = new.u
          a = 1/ (2 * var.eta)
        }
      }
      
      all.var.eta=c(all.var.eta, var.eta)

      if((norm(r.prime,type="2")<eps.prime) & (norm(r.dual,type="2")<eps.dual)){break}

    }
    
    alpha_all[,il]=as.matrix(new.b)

    
    yhat = W%*%new.b 
  
    Hmat=W%*%solve.mat%*%t(W) #
    df=sum(diag(Hmat))

  
    gacv= sum(rhotau((y-yhat),tau))/(length(y)-sum(diag(Hmat))) 

    gacv_all[il] = gacv
    df_all[il]=df
    eta_all[il] = var.eta
  }
  

  j=which.min(gacv_all)
  gacv= gacv_all[j]
  df = df_all[j]
 
  eta=eta_all[j]
  lambdac=lambda[j,]
  alpha_hat=alpha_all[,j]
  gamma=alpha_hat
  
  gamma.mtx=matrix(gamma,J,np)
  theta=Q2%*%gamma.mtx
  beta=B%*%theta
  end_time = Sys.time()
  result_time=end_time - start_time  


  list(gamma = gamma, theta = theta, beta = beta,
       gacv = gacv, 
       lambdac = lambdac, time = result_time, 
       iter = iter.admm)
}

