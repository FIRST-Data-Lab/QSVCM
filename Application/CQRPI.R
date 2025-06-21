source('cv.pred.qsvcm.R')

#' Calculate average coverage probabilities and lengths of 95\% prediction intervals for a new response variable using QBPST-based conformal quantile regression..
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
#' @param lambda The vector of the candidates of penalty parameter.
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
#' @param rate the rate of the training and calibration sets -- default is 0.50.
#' \cr
#' @param cqr.seed seed for conformalized quantile regression -- default is 123.
#' 
#' @return The function returns a list with the following items:
#' \item{coverage.rate}{Average coverage probabilities of 95\% prediction intervals for a new response variable using QBPST-based conformal quantile regression.}
#' \item{mean.coverage.length}{Average lengths of 95\% prediction intervals for a new response variable using QBPST-based conformal quantile regression.}


CQRPI = function(y, X, S, V, Tr, d = 2, r = 1, 
                 lambda = 10^(seq(log10(10^(-6)), log10(10^(1)), length.out = 10)),
                  eta = 1.0, 
                  eps.abs = 10^(-4), eps.rel = 10^(-2), step = 1, mu.eta = 10, incr = 2,
                  B = NULL, Q2 = NULL, P= NULL, 
                  rate = 0.50, cqr.seed = 123)
{

  lambda.cv = lambda
  all.cv.result975 = list()
  all.cv.measure975 = c()
  all.cv.result025 = list()
  all.cv.measure025 = c()
  for(lam_iter in 1:nlambda){
    cvlam = lambda.cv[lam_iter] 
    cvresult975 = cv.pred.qsvcm(y, X, S, V, Tr, d = d, r = r, lambda = cvlam, tau = 0.975, nfold = 5)
    cvresult025 = cv.pred.qsvcm(y, X, S, V, Tr, d = d, r = r, lambda = cvlam, tau = 0.025, nfold = 5)
    
    all.cv.result975[[lam_iter]] = cvresult975
    all.cv.result025[[lam_iter]] = cvresult025
    
    cv.measure975 = c(mean(cvresult975[[1]]), mean(cvresult975[[2]]))
    cv.measure025 = c(mean(cvresult025[[1]]), mean(cvresult025[[2]]))
    
    
    all.cv.measure975 =  rbind(all.cv.measure975, cv.measure975)
    all.cv.measure025 =  rbind(all.cv.measure025, cv.measure025)
  }
  
  
  # set up for smoothing parameters in the penalty term (update)
  initial.cv.result975 = cbind(round(lambda.cv, 2), all.cv.measure975)
  ind.cv975 = which.min(initial.cv.result975[,2])
  initial.lambda975 = lambda.cv[ind.cv975]
  
  initial.cv.result025 = cbind(round(lambda.cv, 2), all.cv.measure025)
  ind.cv025 = which.min(initial.cv.result025[,2])
  initial.lambda025 = lambda.cv[ind.cv025]
  
  new.nlambda = 10 
  lambda.scale = 100
  new.lambda975 = 10^(seq(log10(initial.lambda975/lambda.scale),log10(initial.lambda975*lambda.scale),length.out=new.nlambda))
  new.lambda025 = 10^(seq(log10(initial.lambda025/lambda.scale),log10(initial.lambda025*lambda.scale),length.out=new.nlambda))
  
  
  lambda.cv975 = new.lambda975
  lambda.cv025 = new.lambda025
  
  all.cv.result2.975 = list()
  all.cv.measure2.975 = c()
  all.cv.result2.025 = list()
  all.cv.measure2.025 = c()
  
  for(lam_iter in 1:new.nlambda){
    cvlam975 = lambda.cv975[lam_iter] 
    cvlam025 = lambda.cv025[lam_iter] 
    
    cvresult975 = cv.pred.qsvcm(y, X, S, V, Tr, d = d, r = r, lambda = cvlam975, tau = 0.975, nfold = 5)
    cvresult025 = cv.pred.qsvcm(y, X, S, V, Tr, d = d, r = r, lambda = cvlam025, tau = 0.025, nfold = 5)
    
    all.cv.result2.975[[lam_iter]] = cvresult975
    all.cv.result2.025[[lam_iter]] = cvresult025
    
    cv.measure975 = c(mean(cvresult975[[1]]), mean(cvresult975[[2]]))
    cv.measure025 = c(mean(cvresult025[[1]]), mean(cvresult025[[2]]))
    
    
    all.cv.measure2.975 =  rbind(all.cv.measure2.975, cv.measure975)
    all.cv.measure2.025 =  rbind(all.cv.measure2.025, cv.measure025)
  }
  
  
  final.cv.result975 = cbind(round(lambda.cv975, 2), all.cv.measure2.975)
  final.cv.result025 = cbind(round(lambda.cv025, 2), all.cv.measure2.025)
  
  ind.cv975 = which.min(final.cv.result975[,2])
  ind.cv025 = which.min(final.cv.result025[,2])
  
  
  final.lambda975 = lambda.cv975[ind.cv975]
  final.lambda025 = lambda.cv025[ind.cv025]
  
 
  
  All.result.QRCPI = c()
  All.result.QRCPI.length = c()
  nfold = 5
  n = length(y)
  set.seed(cqr.seed)
  Test = sample(1:n)
  sfold = round(n / nfold)
  
  for(ii in 1:nfold){
    if(ii < nfold){
      Test.set = sort(Test[((ii - 1) * sfold + 1):(ii * sfold)])
    }
    if(ii == nfold){
      Test.set = sort(Test[((ii - 1) * sfold + 1):n])
    }
    Train.set = setdiff(1:n, Test.set)
    
    if(is.vector(X) == 1){
      X.test = as.matrix(X[Test.set])
      X.train = as.matrix(X[Train.set])
    } else {
      X.test = X[Test.set, ]
      X.train = X[Train.set, ]
      S.test = S[Test.set,]
      S.train = S[Train.set, ]
    }
    
    B.test = B[Test.set, ]
    B.train = B[Train.set, ]
    y.test = y[Test.set]
    y.train = y[Train.set]
    
    BQ2.test = as.matrix(B.test%*%Q2) 
    
    
   
    BQ2.new = BQ2.test
    X.new = X.test
    Y.new = y.test

    I1 = sample(1:length(y.train), size = rate * length(y.train))
    I2 = setdiff(1:length(y.train), I1)
    
    y.I1 <- y.train[I1]; X.I1 <- X.train[I1,]; S.I1 <- S.train[I1,]; B.I1 <- B.train[I1,];
    y.I2 <- y.train[I2]; X.I2 <- X.train[I2,]; S.I2 <- S.train[I2,]; B.I2 <- B.train[I2,];
    
    BQ2.I1 <- as.matrix(B.I1%*%Q2) 
    BQ2.I2 <- as.matrix(B.I2%*%Q2) 
    
    W.I1 = as.matrix(kr(X.I1,BQ2.I1,byrow=TRUE))
    W.I2 = as.matrix(kr(X.I2,BQ2.I2,byrow=TRUE))
    np.I1 = ncol(X.I1)
    J.I1 = dim(BQ2.I1)[2]
    
    np.I2 = ncol(X.I2)
    J.I2 = dim(BQ2.I2)[2]
    

    Reference.result025 = fit.qsvcm(y.I1, X.I1, S.I1, V, Tr, d, r, 
                                      final.lambda025, tau = 0.025, eta, eps.abs, eps.rel, step, mu.eta, incr, B.I1, Q2, P)
      
    Reference.result975 = fit.qsvcm(y.I1, X.I1, S.I1, V, Tr, d, r, 
                                      final.lambda975, tau = 0.975, eta, eps.abs, eps.rel, step, mu.eta, incr, B.I1, Q2, P)
    
    
    yhat975.I2 = W.I2%*%Reference.result975$gamma 
    yhat025.I2 = W.I2%*%Reference.result025$gamma 
    
    scores.CQR <- apply(cbind(yhat025.I2-y.I2, y.I2-yhat975.I2), 1, max)
    
    emp.quantile = 0.95*(1+1/length(y.I2)) 
    Q.E.I2 <- quantile(scores.CQR, emp.quantile)
    
    
    
    W.new = as.matrix(kr(X.new, BQ2.new, byrow = TRUE))
    
    yhat975.new = W.new%*%Reference.result975$gamma 
    yhat025.new = W.new%*%Reference.result025$gamma 
    
   
    
    pred.int = cbind(yhat025.new-Q.E.I2, yhat975.new + Q.E.I2)
    Y.new.cover = (pred.int[,1] <= Y.new) & (pred.int[,2] >= Y.new)
    Length.pred = pred.int[,2]-pred.int[,1]
    mean.length = mean(Length.pred)
    cov.rate = mean(Y.new.cover)
    
    
    All.result.QRCPI = c(All.result.QRCPI, cov.rate)
    All.result.QRCPI.length = c(All.result.QRCPI.length, mean.length)
    cat("\nCoverage.Rate:", cov.rate)
    cat("\nCoverage.Mean.Length:", mean.length)
  }
  

  list(coverage.rate = All.result.QRCPI, mean.coverage.length = All.result.QRCPI.length)
}
