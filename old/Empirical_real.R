

library(kernlab)
library(drf)
library(Matrix)
# Load necessary libraries
library(mvtnorm) # for generating multivariate normal random variables
library(ks)
library(doParallel)

library(dplyr)
library(kableExtra)
library(copula)
library("MulvariateRandomForestVarImp")
library(doParallel)


source("drfnew_v2.R")
#source("applications")
source("genData.R")
source("evaluation.R")


cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)






##########################################
### Start Evaluation of wage data #######
##########################################

start.time <- Sys.time()



set.seed(10)
n<-2000
num.trees<-500

res_wage<-evalrealdata(dataset="real_wagedata", n=n, ntest=round(1/2*n), num.trees=num.trees)


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

sort( round(res_wage$ressynth$VI,3), decreasing=T)



# Redefine when data is loaded
set.seed(10)
num.trees<-500
# Redefine when data is loaded

Xtest<-res_wage$Xtest
Ytest<-res_wage$Ytest
X<-res_wage$X
Y<-res_wage$Y
VI<-res_wage$ressynth$VI


metric<-"NPLD"


sort(VI,decreasing=T)

plot(1:length(VI),sort(VI), cex=0.8, col="darkblue", ylab="Sorted VI", xlab="")
#text(x=65, y= VI[ "education_level" ], labels="education level"  )
#text(x=68, y= sort(VI)["age"], labels="age"  )
#text(x=70, y= sort(VI)["male"]+0.03, labels="male"  )



eval<-rep(NA,2)
evalMAD<-rep(NA,2)
p<-ncol(X)

DRFall <- drf(X=X, Y=Y, num.trees=num.trees)

quantpredictall<-predict(DRFall, newdata=Xtest, functional="quantile",quantiles=c(0.5))
evalMAD[1] <- mean(sapply(1:nrow(Xtest), function(j)  abs(Ytest[j] - quantpredictall$quantile[,,"q=0.5"][1]) ))

eval[1]<-distpredicteval(X,Y,Xtest, Ytest,d=metric, num.trees=num.trees, parallel=F)


cutoff<-0.01

DRFred <- drf(X=X[,names(VI[VI > cutoff])], Y=Y, num.trees=num.trees)
quantpredictred<-predict(DRFred, newdata=Xtest[,names(VI[VI > cutoff])], functional="quantile",quantiles=c(0.5))
evalMAD[2] <- mean(sapply(1:nrow(Xtest), function(j)  abs(Ytest[j] - quantpredictred$quantile[,,"q=0.5"][1]) ))

eval[2]<-distpredicteval(X[,names(VI[VI > cutoff])],Y,Xtest[,names(VI[VI > cutoff])], Ytest,d=metric, parallel=F, num.trees=num.trees)

round(sort(VI[VI > cutoff], decreasing=T),3)



print("cutoff=0.01")
round((eval[2]-eval[1])/eval[1]*100,3)
round((evalMAD[2]-evalMAD[1])/evalMAD[1]*100,3)


cutoff<-0.001

DRFred <- drf(X=X[,names(VI[VI > cutoff])], Y=Y, num.trees=num.trees)
quantpredictred<-predict(DRFred, newdata=Xtest[,names(VI[VI > cutoff])], functional="quantile",quantiles=c(0.5))
evalMAD[2] <- mean(sapply(1:nrow(Xtest), function(j)  abs(Ytest[j] - quantpredictred$quantile[,,"q=0.5"][1]) ))

eval[2]<-distpredicteval(X[,names(VI[VI > cutoff])],Y,Xtest[,names(VI[VI > cutoff])], Ytest,d=metric,parallel=F, num.trees=num.trees)

round(sort(VI[VI > cutoff], decreasing=T),3)

print("cutoff=0.001")
round((eval[2]-eval[1])/eval[1]*100,3)
round((evalMAD[2]-evalMAD[1])/evalMAD[1]*100,3)
##########################################
### End Evaluation of wage data #######
##########################################





##########################################
### Start Evaluation of birth data #######
##########################################

start.time<- Sys.time()
set.seed(10)


n<-2000
num.trees<-500

res_birth<-evalrealdata(dataset="real_birthdata", n=n, ntest=round(1/2*n), num.trees=num.trees)


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


Xtest<-res_birth$Xtest
Ytest<-res_birth$Ytest
X<-res_birth$X
Y<-res_birth$Y
VI<-res_birth$ressynth$VI


# Redefine when data is loaded
set.seed(10)
num.trees<-500
# Redefine when data is loaded

metric <- "NPLD"

sort(VI,decreasing=T)

plot(1:length(VI),sort(VI), cex=0.8, col="darkblue", ylab="Sorted VI", xlab="")
#text(x=65, y= VI[ "education_level" ], labels="education level"  )
#text(x=68, y= sort(VI)["age"], labels="age"  )
#text(x=70, y= sort(VI)["male"]+0.03, labels="male"  )


eval<-rep(NA,2)
p<-ncol(X)

DRFall <- drf(X=X, Y=Y, num.trees=num.trees)

eval[1]<-distpredicteval(X,Y,Xtest, Ytest,d=metric, num.trees=num.trees, parallel=T)


cutoff<-0.01

DRFred <- drf(X=X[,names(VI[VI > cutoff])], Y=Y, num.trees=num.trees)
eval[2]<-distpredicteval(X[,names(VI[VI > cutoff])],Y,Xtest[,names(VI[VI > cutoff])], Ytest,d=metric, parallel=T, num.trees=num.trees)


print("cutoff=0.01")
round(sort(VI[VI > cutoff], decreasing=T),3)
round((eval[2]-eval[1])/eval[1]*100,3)



cutoff<-0.001

DRFred <- drf(X=X[,names(VI[VI > cutoff])], Y=Y, num.trees=num.trees)
eval[2]<-distpredicteval(X[,names(VI[VI > cutoff])],Y,Xtest[,names(VI[VI > cutoff])], Ytest,d=metric, parallel=T, num.trees=num.trees)


print("cutoff=0.001")
round(sort(VI[VI > cutoff], decreasing=T),3)
round((eval[2]-eval[1])/eval[1]*100,3)



##########################################
### End Evaluation of birth data #######
##########################################


##########################################
### Start Evaluation of birth data 2 #######
##########################################

start.time <- Sys.time()
set.seed(10)


n<-2000
num.trees<-500

res_birth<-evalrealdata(dataset="real_birthdata2", n=n, ntest=round(1/2*n), num.trees=num.trees)


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


Xtest<-res_birth$Xtest
Ytest<-res_birth$Ytest
X<-res_birth$X
Y<-res_birth$Y
VI<-res_birth$ressynth$VI


# Redefine when data is loaded
set.seed(10)
num.trees<-500
# Redefine when data is loaded

metric <- "MMD"

sort(VI,decreasing=T)

plot(1:length(VI),sort(VI), cex=0.8, col="darkblue", ylab="Sorted VI", xlab="")
#text(x=65, y= VI[ "education_level" ], labels="education level"  )
#text(x=68, y= sort(VI)["age"], labels="age"  )
#text(x=70, y= sort(VI)["male"]+0.03, labels="male"  )


eval<-rep(NA,2)
p<-ncol(X)

DRFall <- drf(X=X, Y=Y, num.trees=num.trees)

eval[1]<-distpredicteval(X,Y,Xtest, Ytest,d=metric, num.trees=num.trees, parallel=T)


cutoff<-0.01

DRFred <- drf(X=X[,names(VI[VI > cutoff])], Y=Y, num.trees=num.trees)
eval[2]<-distpredicteval(X[,names(VI[VI > cutoff])],Y,Xtest[,names(VI[VI > cutoff])], Ytest,d=metric, parallel=T, num.trees=num.trees)


print("cutoff=0.01")
round(sort(VI[VI > cutoff], decreasing=T),3)
round((eval[2]-eval[1])/eval[1]*100,3)



cutoff<-0.001

DRFred <- drf(X=X[,names(VI[VI > cutoff])], Y=Y, num.trees=num.trees)
eval[2]<-distpredicteval(X[,names(VI[VI > cutoff])],Y,Xtest[,names(VI[VI > cutoff])], Ytest,d=metric, parallel=T, num.trees=num.trees)


print("cutoff=0.001")
round(sort(VI[VI > cutoff], decreasing=T),3)
round((eval[2]-eval[1])/eval[1]*100,3)



##########################################
### End Evaluation of birth data 2 #######
##########################################










### Testing of Testing

n<-2000
tmp<-genData(dataset = "bivariatesynthetic", n = n, p = 4)


X<-tmp$X
Y<-as.matrix(tmp$y)
colnames(X) <- paste0("X",1:ncol(X))

B<-100
num.trees=500
num.features=10
sample.splitting=T
ntest=10#round(0.2*n)


ressynth<-drfwithVI(X, Y, B=B, num.trees=num.trees, num.features=num.features, sample.splitting=sample.splitting, ntest=ntest)


ressynth$VI
ressynth$VIcorrected

ressynth$I0list[[1]]
ressynth$I0list[[2]]




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
