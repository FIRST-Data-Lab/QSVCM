source('QSVCM_est.R') 
#' Fitting quantile spatially varying coefficient models
#'
#' \code{fit.qsvcm} fits the quantile spatially varying coeffcient models. See \code{QSVCM_est} for detailed arguments.
#'
fit.qsvcm =
  function(y, X, S, V, Tr, d = 2, r = 1, lambda = 10^(seq(log10(10^(-6)), log10(10^(1)), length.out = 10)), tau = 0.5, eta = 1.0, 
           eps.abs = 10^(-4), eps.rel = 10^(-2), step = 1, mu.eta = 10, incr = 2,
           B = NULL, Q2 = NULL, P= NULL)
{
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
    
    this.call = match.call()

    if ((is.null(B))|(is.null(Q2))|(is.null(P))){
      Basis.full = basis(V, Tr, d, r, S)
      K = Basis.full$K
      Q2 = Basis.full$Q2
      B = Basis.full$B
      ind.inside = Basis.full$Ind.inside
      tria.all = Basis.full$tria.all
    
      y = y[ind.inside]
      X = as.matrix(X[ind.inside, ])
      S = as.matrix(S[ind.inside, ])
    
      BQ2 = B %*% Q2
      P = t(Q2) %*% K %*% Q2
    }
    mfit = QSVCM_est(y = y, X = X, B = B, Q2 = Q2, P = P, 
                     lambda = lambda, 
                     tau = tau, eta = eta, eps.abs = eps.abs, eps.rel = eps.rel,
                     step = step, mu.eta = mu.eta, incr = incr)
    
    mfit$V = V; mfit$Tr = Tr; mfit$d = d; mfit$r = r; mfit$B = B;
    mfit$Q2 = Q2; mfit$P = P; mfit$K = K;  mfit$X = X; mfit$y = y; mfit$S = S;
    if ((is.null(B))|(is.null(Q2))|(is.null(P))){
           mfit$ind.inside = ind.inside;
    mfit$tria.all = tria.all;
     }
    mfit$call = this.call;
    class(mfit) = "qsvcm"
    return(mfit)
    
}