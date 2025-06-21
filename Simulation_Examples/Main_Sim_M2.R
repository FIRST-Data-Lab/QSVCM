rm(list = ls())

# install.packages("devtools",repos="http://cran.r-project.org")
# library(devtools)
# install_github("funstatpackages/BPST")

#library
library('mgcv')
library('MGLM') 
library('BPST')

#source functions
source('stoper.R')
source('QSVCM_est.R') 
source('rhotau.R')
source('QmixNorm.R')
source('QRWboot.R')  
source('fit.qsvcm.R')  
source('CQRPI.R')

n = 1000     # number of sample points
tau = 0.50  # quantile levels
dist = "MN" ## choose "MN" for normal mixture distribution 
            ## choose "t2" for t2 error distribution 

if (dist == "MN"){
  Finv = QmixNorm(q = tau, mu = c(-5, 0, 5), sd = sqrt(c(4, 0.5, 4)), pi = c(0.1, 0.8, 0.1))
  pop.all=as.matrix(read.csv('eg2_pop_mnorm.csv', header=TRUE))    
} else if(dist == "t2"){
  Finv = qt(tau,2)
  pop.all=as.matrix(read.csv('eg2_pop_t2.csv', header=TRUE))
} else{
  stop("Choose MN or t2 for error distribution.")
}


## rough triangulation (Tr 1)
Tr=as.matrix(read.csv('Tr_1.csv',header=FALSE))
V=as.matrix(read.csv('V_1.csv',header=FALSE))

## fine triangulation (Tr 2)
# Tr=as.matrix(read.csv('Tr_2.csv',header=FALSE)) 
# V=as.matrix(read.csv('V_2.csv',header=FALSE))

run.bootstrap = 0    # 1: run bootstrap procedure for confidence interval 
# 0: skip bootstrap results

run.CQRPI = 0 # 1: run CQR prediction intervals
# 0: skip CPR prediciton intervals

nsim=100


 
# Population
N.all=nrow(pop.all)
ind1=inVT(V,Tr,pop.all[,'u'],pop.all[,'v'])
ind1=ind1$ind.inside
ind2=(1:nrow(pop.all))[!is.na(pop.all[,1])]
ind=sort(intersect(ind1,ind2))
pop.r=pop.all[ind,]
Y.pop=pop.r[,1]


beta.pop=pop.r[,c('m1','m2','m3','m4')]
beta.pop[,1] = beta.pop[,1] + Finv*beta.pop[,3]
beta.pop[,2] = beta.pop[,2] + Finv*beta.pop[,4]


X.pop=pop.r[,c('x1','x2')]
S.pop=round(pop.r[,c('u','v')],2)
Npop=nrow(pop.r)
np=ncol(X.pop)

u=unique(round(pop.all[,'u'],2))
v=unique(round(pop.all[,'v'],2))

n1=length(u)
n2=length(v)

# Population basis
d=3; r=1;
B0.pop=basis(V,Tr,d,r,S.pop)
Q2=B0.pop$Q2
B.pop=B0.pop$B
BQ2.pop=as.matrix(B.pop%*%Q2)
K=B0.pop$K
P=t(Q2)%*%K%*%Q2



lambda_start=log10(10^(-6))
lambda_end=log10(10^(1))
nlambda=10
lambda=10^(seq(lambda_start,lambda_end,length.out=nlambda))



J=dim(Q2)[2]

MISE1_all=c()
MISE2_all=c()
mhat_all_all=c()

all_all_Beta1 = c()
all_all_Beta2 = c()

all_all_cp = c()
all_all_length = c()


elapsed.time_all=c()
for(iter in 1:nsim){
  cat("\nIteration:",iter)
  set.seed(iter)
  
  ind.s=sort(sample(Npop,n,replace=FALSE))
  data=as.matrix(pop.r[ind.s,])
  Y=data[,1]
  beta0=data[,c('m1','m2','m3','m4')]
  X=data[,c('x1','x2')]
  S=data[,c('u','v')]
  
  # generate Bivariate spline basis 
  B=B.pop[ind.s,]
  
  
  # sample population basis. 
  true.alpha.samp = beta0
  true.alpha.samp[,1] = true.alpha.samp[,1] + Finv*true.alpha.samp[,3]
  true.alpha.samp[,2] = true.alpha.samp[,2] + Finv*true.alpha.samp[,4]
 
  
  new.result=QSVCM_est(Y, X, B, Q2, P, lambda, tau)
  initial.lambda=new.result$lambdac
  
  new.nlambda=10
  lambda.scale=100
  new.lambda1 = 10^(seq(log10(initial.lambda[1]/lambda.scale),log10(initial.lambda[1]*lambda.scale),length.out=nlambda))
  new.lambda2 = 10^(seq(log10(initial.lambda[2]/lambda.scale),log10(initial.lambda[2]*lambda.scale),length.out=nlambda))
  new.lambda = cbind(new.lambda1, new.lambda2)

  new.result2 = QSVCM_est(Y, X, B, Q2, P, new.lambda, tau)
  new.b=new.result2$gamma
  final.lambda=new.result2$lambdac
  
  
  beta_hat=new.result2$beta
  beta.r1=BQ2.pop%*%new.b[1:J]
  beta.r2=BQ2.pop%*%new.b[(J+1):(2*J)]

  mhat_all=cbind(beta.r1,beta.r2)

  
  # Confidence intervals for population grids
  if (run.bootstrap == 1){
    result_CB = QRWboot(y = Y, X = X, S = S, V = V, Tr = Tr, d = d, r = r, 
                        BQ2_pop = BQ2.pop, tau = tau, nB = 100)
    
    all_Beta1 = (result_CB$Lower.CB[,1]<beta.pop[,1])&(result_CB$Upper.CB[,1]>beta.pop[,1])
    all_Beta2 = (result_CB$Lower.CB[,2]<beta.pop[,2])&(result_CB$Upper.CB[,2]>beta.pop[,2])
    all_all_Beta1 = cbind(all_all_Beta1, all_Beta1)
    all_all_Beta2 = cbind(all_all_Beta2, all_Beta2)
  } 
  
  # CQR for a new response
  if (run.CQRPI == 1){
    result_CQRPI = CQRPI(y = Y, X = X, S = S, V = V, Tr = Tr, d = d, r = r, B = B, Q2 = Q2, P= P)
    
    coverage.rate = result_CQRPI$coverage.rate
    mean.coverage.length = result_CQRPI$mean.coverage.length

	all_all_cp = c(all_all_cp,  mean(coverage.rate))
	all_all_length = c(all_all_length, mean(mean.coverage.length))

  } 
  
 
  # MISE
  MISE1=(0.02*0.02)*sum((beta.pop[,1]-mhat_all[,1])^2, na.rm=TRUE)
  MISE2=(0.02*0.02)*sum((beta.pop[,2]-mhat_all[,2])^2, na.rm=TRUE)
  
 
  MISE1_all=cbind(MISE1_all,MISE1)
  MISE2_all=cbind(MISE2_all,MISE2)

  
  mhat_all_all=mhat_all 
  
  cat("\nMISE1:", MISE1,"MISE2:", MISE2)
}  

 
RMISE1_all=sqrt(MISE1_all)
RMISE2_all=sqrt(MISE2_all)
 
RMISE1=mean(RMISE1_all)
RMISE2=mean(RMISE2_all)

SD_RMISE1 = sd(RMISE1_all)
SD_RMISE2 = sd(RMISE2_all)

result2=cbind(RMISE1, RMISE2, SD_RMISE1, SD_RMISE2)
result2
 



if (run.bootstrap == 1){
  cat("95% Coverage Percentage(Beta1):", mean(apply(all_all_Beta1, 1, mean)), "\n");  cat("95% Coverage Percentage(Beta2):", mean(apply(all_all_Beta2, 1, mean)))
}


if (run.CQRPI == 1){
  cat("Coverage probaility:", mean(all_all_cp), "\n");  cat("Average length:", mean(all_all_length))
}



# ######################################################################
# # Only for the plots
# ######################################################################
# library(plot3D)
# # true surface
# beta.r1.all=matrix(NA,n1*n2,1)
# beta.r1.all[ind]= as.matrix(beta.pop[,1])
# beta.r1.mtx=matrix(beta.r1.all,n1,n2)
# beta.r2.all=matrix(NA,n1*n2,1)
# beta.r2.all[ind]= as.matrix(beta.pop[,2])
# beta.r2.mtx=matrix(beta.r2.all,n1,n2)
# 
# # true surfaces of beta0 and beta1
# image2D(u,v,z=beta.r1.mtx, col = heat.colors(5))
# contour(u,v,matrix(beta.r1.mtx,n1,n2), add = TRUE, drawlabels = TRUE)
# image2D(u,v,z=beta.r2.mtx, col = heat.colors(5))
# contour(u,v,matrix(beta.r2.mtx,n1,n2), add = TRUE, drawlabels = TRUE)
#  
# ## estimated surfaces of beta0 and beta1
# beta.r1.all=matrix(NA,n1*n2,1)
# beta.r1.all[ind,]=mhat_all[,1]
# beta.r1.mtx=matrix(beta.r1.all,n1,n2)
# image2D(u,v,z=matrix(beta.r1.mtx,n1,n2), col = heat.colors(5))
# contour(u,v,beta.r1.mtx, add = TRUE, drawlabels = TRUE)
# 
# beta.r2.all=matrix(NA,n1*n2,1)
# beta.r2.all[ind,]=mhat_all[,2]
# beta.r2.mtx=matrix(beta.r2.all,n1,n2)
# image2D(u,v,z=matrix(beta.r2.mtx,n1,n2), col = heat.colors(5))
# contour(u,v,beta.r2.mtx, add = TRUE, drawlabels = TRUE)