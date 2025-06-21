##############################################################################
# Main_Application.R
##############################################################################
rm(list = ls())

#library
# install.packages("devtools",repos="http://cran.r-project.org")
# library(devtools)
# install_github("funstatpackages/BPST")

library('BPST')
library('mgcv') 
library('MGLM') 

#source functions
source('stoper.R') 
source('fit.qsvcm.R')
source('rhotau.R')
source("QRWboot.R")
source('cv.pred.qsvcm.R')
source('CQRPI.R')

# Quantile level
tau = 0.50 

# degree and smoothness
d = 4; r = 1

# Load matrices for triangles (Tr) and vertices (V)
Tr=as.matrix(read.csv('Tr_usa.csv', header = F))
V=as.matrix(read.csv('V_usa.csv', header = F))

# Load data
Mort_data = read.csv('Data_mortality.csv', header = T)
# Generating population grids
S_pop = as.matrix(read.csv('pop_location_d025.csv', header = T)) 

# Population basis and index
B0_pop = basis(V,Tr,d,r,S_pop)
Q2 = B0_pop$Q2
nq = dim(Q2)[2]
BQ2_pop = as.matrix(B0_pop$B%*%Q2)
Ind_all = B0_pop$Ind.inside
K = B0_pop$K
P = t(Q2)%*%K%*%Q2


# Load response and covariates / standardized covariates
Y = as.matrix(Mort_data[, 7])
X = Mort_data[, 1:4]
X = as.matrix(scale(X, scale = TRUE)) 
X = as.matrix(cbind(1, X))
S = cbind(Mort_data$Longtitude, Mort_data$Latitude)

# Generate sample Bivariate spline basis 
B0 = basis(V, Tr, d, r, S)
Ind = B0$Ind.inside
Y = Y[Ind]; X = X[Ind,]; S = S[Ind,]
B = B0$B
BQ2 = as.matrix(B%*%Q2) 


# set up for smoothing parameters in the penalty term
lambda_start = log10(10^(-6))
lambda_end = log10(10^(1))
nlambda = 10
lambda = 10^(seq(lambda_start, lambda_end, length.out = nlambda))


lambda.cv = lambda
all.cv.result = list()
all.cv.measure = c()
for(lam_iter in 1:nlambda){
  cvlam = lambda.cv[lam_iter] 
  cvresult = cv.pred.qsvcm(Y, X, S, V, Tr, d = d, r = r, lambda = cvlam, tau = tau, nfold = 5)
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
  cvresult = cv.pred.qsvcm(Y, X, S, V, Tr, d = d, r = r, lambda = cvlam, tau = tau, nfold = 5)
  all.cv.result2[[lam_iter]] = cvresult
  cv.measure = c(mean(cvresult[[1]]), mean(cvresult[[2]]))
  all.cv.measure2 =  rbind(all.cv.measure2, cv.measure)
}


final.cv.result = cbind(round(lambda.cv, 2), all.cv.measure2)
ind.cv = which.min(final.cv.result[,2])
final.lambda = lambda.cv[ind.cv]

# 5-cv results
final.cv.result[ind.cv,]

# fit the final model
new.result2 = fit.qsvcm(y = Y, X = X, S = S, V = V, Tr =Tr, d = d, r = r,
                        lambda = final.lambda, tau = tau)
result = new.result2
y_hat = kr(X, BQ2, byrow = TRUE)%*%result$gamma

# obtain estimates for bivariate functions on the domain:
np=ncol(X)
mhat=matrix(0,dim(BQ2)[1],np)
Jn=dim(BQ2)[2]
mhat=c()
mhat_all=c()
for(i in 1:np){
  mhat = cbind(mhat,BQ2%*%result$gamma[(1+(i-1)*Jn):(Jn*i)])
  mhat_all = cbind(mhat_all,BQ2_pop%*%result$gamma[(1+(i-1)*Jn):(Jn*i)])
}

 
# MQLPE & MSPE
mean(rhotau((Y-y_hat),tau=tau)); mean((Y-y_hat)^2)
 


# Confidence intervals for population grids
run.bootstrap = 0     # 1: run bootstrap procedure for confidence interval 
           # 0: skip bootstrap results
if (run.bootstrap == 1){
  result_CB = QRWboot(y = Y, X = X, S = S, V = V, Tr = Tr, d = d, r = r, BQ2_pop = BQ2_pop, tau = tau, nB = 100)

  Lower_mhat_all = result_CB$Lower.CB
  Upper_mhat_all = result_CB$Upper.CB
} 


all_all_cp  = c()
all_all_length = c()

# CQR prediction intervals for a new response:
run.CQRPI = 0 # 1: run CQR prediction intervals
# 0: skip CPR prediction intervals

  # CQR for a new response
  if (run.CQRPI == 1){
    result_CQRPI = CQRPI(y = Y, X = X, S = S, V = V, Tr = Tr, d = d, r = r, B = B, Q2 = Q2, P= P)
    
    coverage.rate = result_CQRPI$coverage.rate
    mean.coverage.length = result_CQRPI$mean.coverage.length

   	all_all_cp = c(all_all_cp,  mean(coverage.rate))
  	all_all_length = c(all_all_length, mean(mean.coverage.length))

  } 


 