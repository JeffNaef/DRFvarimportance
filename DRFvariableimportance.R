

library(kernlab)
library(drf)
library(Matrix)
# Load necessary libraries
library(mvtnorm) # for generating multivariate normal random variables
library(ks)

source("drfnew_v2.R")
#source("applications")
source("genData.R")

set.seed(10)

### Continue with more fancy examples!!!

n<-200
ntest<-round(n*0.1)
num.trees<-100

## Step 1: Get the dataset
tmp<-genData(dataset = "synthetic1", n = n, p = 10, meanShift = 1, sdShift = 1)

X<-tmp$X
Y<-as.matrix(tmp$y)
colnames(X) <- paste0("X",1:ncol(X))




## Step 2: Do Analysis

### 2.a) if the dataset is synthetic, we can check the correct variable ordering
ressynth<-drfwithVI(X, Y, B=1, num.trees=num.trees)

sort(ressynth$VI)


### 2.b) if it is a real dataset, I need to make a prediction function that successively filters
# Sample Splitting
Xtest <- X[(round(n - ntest) + 1):n, , drop = F]
Ytest <- Y[(round(n - ntest) + 1):n, , drop = F]
#
X <- X[1:round(n - ntest), , drop = F]
Y <- Y[1:round(n - ntest), , drop = F]



whichvar<-rep(NA,ncol(X)-1)
eval<-rep(NA,ncol(X)-1)
p<-ncol(X)

for (j in 1:(p - 1)  ){
  

  # remove variable with smallest VI
  resreal<-featureeliminnation(X,Y)
  X<-resreal$Xnew
  Xtest<-Xtest[,!( colnames(Xtest) %in% resreal$which), drop=F]
  

  whichvar[j] <-resreal$which
  eval[j]<-distpredicteval(X,Y,Xtest, Ytest,dNPLD, num.trees=num.trees)
}




# # Generate random data for independent variables X
# n <- 1000 # number of observations
# X <-
#   matrix(runif(n * 3), ncol = 3) # 100 x 3 matrix with random values from a standard normal distribution
# 
# d <- 1
# # There is a simple and a fancy example
# example <- "fancy"
# # Sample splitting or not?
# sample.splitting <- T
# ntest <- n/2
# 
# 
# # Define coefficients for the linear combination
# # two-dimensional
# 
# if (example == "simple") {
#   if (d == 2) {
#     B <-
#       matrix(c(2,-1, 0.5, 0, 0, 0), ncol = 3) # 2 x 3 matrix with chosen coefficients
#     mu <- c(0, 0) # mean for the random noise
#     Sigma <-
#       matrix(c(1, 0.5, 0.5, 1), nrow = 2) # covariance matrix for the random noise
#     epsilon <-
#       rmvnorm(n, mu, Sigma) # generate random noise following a bivariate normal distribution
#     
#     # Create the 2-dimensional response variable Y
#     Y <-
#       as.matrix(X %*% t(B) + epsilon) # compute Y as the linear combination of X and B plus random noise
#     
#     
#     
#   } else if (d == 1) {
#     # one-dimensional
#     b <- matrix(c(6, 3, 0), ncol = 3)
#     epsilon <- rnorm(n)
#     ## Create the 2-dimensional response variable Y
#     Y <-
#       as.matrix(X %*% t(b) + epsilon) # compute Y as the linear combination of X and B plus random noise
#     
#     
#     
#   }
# } else if (example == "fancy") {
#   # more fancy example
#   sigX = X[, 2]
#   Y <- matrix(rnorm(n, mean = 0, sd = sqrt(sigX)), nrow = n)
# }

# 
# 
# 
# if (sample.splitting == T) {
#   # Sample Splitting
#   Xtest <- X[(round(n - ntest) + 1):n, , drop = F]
#   Ytest <- Y[(round(n - ntest) + 1):n, , drop = F]
#   #
#   X <- X[1:round(n - ntest), , drop = F]
#   Y <- Y[1:round(n - ntest), , drop = F]
# } else{
#   # No sample splitting
#   Xtest <- X
#   Ytest <- Y
# }
# 
# 
# 
# 
# 
# bandwidth_Y <- drf:::medianHeuristic(Ytest)
# k_Y <- rbfdot(sigma = bandwidth_Y)
# K <- kernelMatrix(k_Y, Y, y = Y)
# 
# 
# if (B > 1){
# 
# DRF <-
#   drfCI(
#     X = X,
#     Y = Y,
#     B = B,
#     num.trees = num.trees
#   )
# 
# # Prediction with all X
# DRFpred = predictdrf(DRF, x = Xtest)
# wall <- DRFpred$weights
# ##
# 
# 
# I0 <- sapply(1:ncol(Xtest), function(j) {
#   # iterate over class 1
#   
#   ## With CI
#   DRFj <-
#     drfCI(
#       X = X[, -j],
#       Y = Y,
#       B = B,
#       num.trees = num.trees
#     )
#   
#   DRFpredj = predictdrf(DRFj, x = Xtest[, -j])
#   wj <- DRFpredj$weights
#   val <- sum(diag( (wj - wall) %*% K %*% t(wj - wall) ))
#   
#   
# 
#     # Get null distribution if B > 1
#     nulldist <- sapply(1:B, function(b) {
#       # iterate over class 1
#       
#       wbj <- DRFpredj$weightsb[[b]]
#       wb <- DRFpred$weightsb[[b]]
#       
#       sum(diag((wb - wall - (wbj - wj)) %*% K %*% t(wb - wall - (wbj - wj))))
#     })
#     ##
#     right_quantile <- quantile(nulldist, 1 - alpha)
#     
#     
#     max(val - unname(right_quantile), 0)
# 
#   
#   
# })
# 
# }else{
#   
#   DRF2<-drf(X,Y, num.trees=num.trees)
#   wall<-predict(DRF2,x = Xtest)$weights
#   
#   I0 <- sapply(1:ncol(Xtest), function(j) {
#     # iterate over class 1
#     
#     ## With CI
#     DRFj <-
#       drf(
#         X = X[, -j],
#         Y = Y,
#         num.trees = num.trees
#       )
#     
#     DRFpredj = predict(DRFj, x = Xtest[, -j])
#     wj <- DRFpredj$weights
#     val <- sum(diag( (wj - wall) %*% K %*% t(wj - wall) ))
#     
#     })
#   
# 
# }
# 
# ### Calculate average weight
# wbar <- colMeans(wall)
# # alternative:
# #wbar<- rep(1/ntest,ntest)
# wall_wbar<-sweep(wall, 2, wbar, "-")
# #( I<-I0/as.numeric( mean(diag(  as.matrix(wall) %*% K %*% t( as.matrix(wall)) )) - colMeans(wall)%*%K%*%colMeans(wall) ) )
# (I <-
#     I0 / as.numeric(sum(diag(
#       wall_wbar %*% K %*% t(wall_wbar)
#     ))))
