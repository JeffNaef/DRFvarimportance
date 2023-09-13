

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
library(kernlab)


source("drfnew_v2.R")
#source("applications")
source("genData.R")
source("evaluation.R")
source("compute_drf_vimp.R")


cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)






##########################################
### Start Evaluation of wage data #######
##########################################

start.time <- Sys.time()



set.seed(10)
n<-2000
ntest <- round(1/2*n)
num.trees<-500

set.seed(10)

tmp<-genData(dataset = "real_wagedata", n = 2*n)

#tmp<-genData(dataset = "motivatingexample", n = n)


X<-tmp$X
Y<-as.matrix(tmp$y)

if (is.null(colnames(X))){
  colnames(X) <- paste0("X",1:ncol(X))
}


Xtest <- X[(round(n - ntest) + 1):n, , drop = F]
Ytest <- Y[(round(n - ntest) + 1):n, , drop = F]
#
X <- X[1:round(n - ntest), , drop = F]
Y <- Y[1:round(n - ntest), , drop = F]


VIMPDRF<-compute_drf_vimp(X, Y, X_test = NULL, num.trees = 500, silent = FALSE)
(VIMPDRF<-sort(VIMPDRF))

resDRF<-evalfinal(VIMPDRF, X ,Y ,Xtest, Ytest, metrics=c("MMD","NPLD","MAD"), num.trees=500 )



# vimp drf
DRF <- drf(X, Y, num.trees = 500)
vimp_drf <- variable_importance(DRF)
vimp_drf<-sort(c(vimp_drf))
names(vimp_drf)<- colnames(X)

resDRF_native<-evalfinal(vimp_drf, X ,Y ,Xtest, Ytest, metrics=c("MMD","NPLD","MAD"), num.trees=500 )




# Sobol MDA
sobolMDA <- list()
for (j in 1:ncol(Y)){
  XY <- as.data.frame(cbind(X, Y[, j]))
  colnames(XY) <- c(paste('X', 1:(ncol(XY)-1), sep=''), 'Y')
  num.trees <- 500
  clock <- Sys.time()
  forest <- sobolMDA::ranger(Y ~., data = XY, num.trees = num.trees, importance = 'sobolMDA')
  duration <- Sys.time() - clock
  print(duration)
  forest$r.squared
  sobolMDA[[j]] <- forest$variable.importance
}

sobolMDA<-colMeans(do.call(rbind,sobolMDA))


ressobolMDA<-evalfinal(sobolMDA, X ,Y ,Xtest, Ytest, metrics=c("MMD","NPLD","MAD"), num.trees=500 )


save(resDRF, resDRF_native, ressobolMDA, n, ntest, X,Y,Xtest,Ytest, file=paste0("real_wagedata_n=", n))



plot(1:ncol(X),resDRF$evalNPLD, type="l", cex=0.5)
plot(1:ncol(X),resDRF$evalMMD, type="l", cex=0.5)
plot(1:ncol(X),resDRF$evalMMD, type="l", cex=0.5)


plot(1:ncol(X),resDRF_native$evalNPLD, type="l", cex=0.5)
plot(1:ncol(X),resDRF_native$evalMMD, type="l", cex=0.5)



plot(1:ncol(X),ressobolMDA$evalNPLD, type="l", cex=0.5)
plot(1:ncol(X),ressobolMDA$evalMMD, type="l", cex=0.5)





##########################################
### Start Evaluation of birth data #######
##########################################




set.seed(10)
n<-2000
ntest <- round(1/2*n)
num.trees<-500

set.seed(10)

tmp<-genData(dataset = "real_birthdata2", n = 2*n)

#tmp<-genData(dataset = "motivatingexample", n = n)


X<-tmp$X
Y<-as.matrix(tmp$y)

if (is.null(colnames(X))){
  colnames(X) <- paste0("X",1:ncol(X))
}


Xtest <- X[(round(n - ntest) + 1):n, , drop = F]
Ytest <- Y[(round(n - ntest) + 1):n, , drop = F]
#
X <- X[1:round(n - ntest), , drop = F]
Y <- Y[1:round(n - ntest), , drop = F]


VIMPDRF<-compute_drf_vimp(X, Y, X_test = NULL, num.trees = 500, silent = FALSE)
(VIMPDRF<-sort(VIMPDRF))

resDRF<-evalfinal(VIMPDRF, X ,Y ,Xtest, Ytest, metrics=c("MMD"), num.trees=500 )



# vimp drf
DRF <- drf(X, Y, num.trees = 500)
vimp_drf <- variable_importance(DRF)
vimp_drf<-sort(c(vimp_drf))
names(vimp_drf)<- colnames(X)

resDRF_native<-evalfinal(vimp_drf, X ,Y ,Xtest, Ytest, metrics=c("MMD"), num.trees=500 )




# Sobol MDA
sobolMDA <- list()
for (j in 1:ncol(Y)){
  XY <- as.data.frame(cbind(X, Y[, j]))
  colnames(XY) <- c(paste('X', 1:(ncol(XY)-1), sep=''), 'Y')
  num.trees <- 500
  clock <- Sys.time()
  forest <- sobolMDA::ranger(Y ~., data = XY, num.trees = num.trees, importance = 'sobolMDA')
  duration <- Sys.time() - clock
  print(duration)
  forest$r.squared
  sobolMDA[[j]] <- forest$variable.importance
}

sobolMDA<-colMeans(do.call(rbind,sobolMDA))


ressobolMDA<-evalfinal(sobolMDA, X ,Y ,Xtest, Ytest, metrics=c("MMD"), num.trees=500 )


save(resDRF, resDRF_native, ressobolMDA, n, ntest, X,Y,Xtest,Ytest, file=paste0("real_birthdata2_n=", n))



#plot(1:ncol(X),resDRF$evalNPLD, type="l", cex=0.5)
plot(1:ncol(X),resDRF$evalMMD, type="l", cex=0.5)
#plot(1:ncol(X),resDRF$evalMAD, type="l", cex=0.5)


#plot(1:ncol(X),resDRF_native$evalNPLD, type="l", cex=0.5)
plot(1:ncol(X),resDRF_native$evalMMD, type="l", cex=0.5)



#plot(1:ncol(X),ressobolMDA$evalNPLD, type="l", cex=0.5)
plot(1:ncol(X),ressobolMDA$evalMMD, type="l", cex=0.5)

















