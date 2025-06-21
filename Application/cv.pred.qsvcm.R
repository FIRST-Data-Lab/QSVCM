#' k-fold cross-validation MQPLPE and MSPE for spatially varying coefficient regression using the BPST
#'
#' This function implements k-fold cross-validation MSPE for spartially varying coefficient regression, and returns the mean squared prediction error.
#'
#' @importFrom BPST basis
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
#' @param lambda The vector of the candidates of penalty parameter -- default is grid points of 10 to the power of a sequence from -6 to 6 by 0.5.
#' \cr
#' @param nfolds The number of folds -- default is 10. Although \code{nfold} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable for \code{nfolds} is 3.
#' \cr
#' @param initial The seed used for cross-validation -- default is 123.
#' \cr
#' @param tau quantile level -- default is 0.50.
#' \cr
#' @return The mean quantile loss prediction error (MQLPE) and the mean squared prediction error (MSPE) based on k-fold cross-validation
#' @details
#'
#'
#' @export
#'
cv.pred.qsvcm =
function(y, X, S, V, Tr, d = 2, r = 1, lambda = 10^(seq(log10(10^(-6)), log10(10^(1)), length.out = 10)), 
         nfold = 10, initial = 123, tau = 0.50
         )
{
  if(nfold < 3){
    warning("The number of folds in CV is too small. Instead, the default 10-fold CV is used.")
    nfold = 10
  }
  if(!is.vector(y)){
    warning("The response variable, y, should be a vector.")
    y = as.vector(y)
  }
  if(!is.matrix(X)){
    warning("The explanatory variable, X, should be a matrix.")
    X = as.matrix(X)
  }
  if(!is.matrix(S)){
    warning("The coordinates, S, should be a matrix.")
    S = as.matrix(S)
  }

  Ball = basis(V, Tr, d, r, S)
  K = Ball$K
  Q2 = Ball$Q2
  B = Ball$B
  ind.inside = Ball$Ind.inside
  tria.all = Ball$tria.all
  
  ################ ind.inside sample #################

  P=t(Q2)%*%K%*%Q2
  y = y[ind.inside]
  X = X[ind.inside, ]
  S = S[ind.inside, ]
  
  n = length(y)
  sfold = round(n / nfold)
  set.seed(initial)
  Test = sample(1:n)
  cv.error1 = c()
  cv.error2 = c()
  cv.error3 = c()
  cv.error4 = c()

  cv.coverp = c()
  cv.lengthp = c()
  all.lambda = c()

  
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
    }
    
    
    B.test = B[Test.set, ]
    B.train = B[Train.set, ]
    y.test = y[Test.set]
    y.train = y[Train.set]
    
    mfit.ii = QSVCM_est(y.train, X.train, B.train, Q2, P, lambda, tau = tau)
    W.test = as.matrix(kr(X.test, B.test%*%Q2, byrow = TRUE))  
    ypred.ii = W.test %*%  as.vector(mfit.ii$gamma)

    
    ## quantile regression.
    pred.error1 = mean(rhotau(y.test - ypred.ii, tau = tau)) 
    pred.error2 = mean((y.test - ypred.ii)^2)
   

    cv.error1 = c(cv.error1, pred.error1)
    cv.error2 = c(cv.error2, pred.error2)

  }
  list(cv.error1, cv.error2) 
}
