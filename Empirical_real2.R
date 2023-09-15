

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

B<-10

resDRF<-list()
resDRF_native<-list()
ressobolMDA<-list()

tmp<-genData(dataset = "real_wagedata", n = n)

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



for (b in 1:B){




VIMPDRF<-compute_drf_vimp(X, Y, X_test = NULL, num.trees = 500, silent = FALSE)
(VIMPDRF<-sort(VIMPDRF))

resDRF[[b]]<-evalfinal(VIMPDRF, X ,Y ,Xtest, Ytest, metrics=c("MMD","NPLD","MAD"), num.trees=500 )



# vimp drf
DRF <- drf(X, Y, num.trees = 500)
vimp_drf <- variable_importance(DRF)
vimp_drf<-sort(c(vimp_drf))
names(vimp_drf)<- colnames(X)

resDRF_native[[b]]<-evalfinal(vimp_drf, X ,Y ,Xtest, Ytest, metrics=c("MMD","NPLD","MAD"), num.trees=500 )




# Sobol MDA
sobolMDAj <- list()
for (j in 1:ncol(Y)){
  XY <- as.data.frame(cbind(X, Y[, j]))
  colnames(XY) <- c(paste('X', 1:(ncol(XY)-1), sep=''), 'Y')
  num.trees <- 500
  clock <- Sys.time()
  forest <- sobolMDA::ranger(Y ~., data = XY, num.trees = num.trees, importance = 'sobolMDA')
  duration <- Sys.time() - clock
  print(duration)
  forest$r.squared
  sobolMDAj[[j]] <- forest$variable.importance
}

sobolMDA<-colMeans(do.call(rbind,sobolMDAj))


ressobolMDA[[b]]<-evalfinal(sobolMDA, X ,Y ,Xtest, Ytest, metrics=c("MMD","NPLD","MAD"), num.trees=500 )
}

save(resDRF, resDRF_native, ressobolMDA, X,Y, Xtest, Ytest, n, ntest, file=paste0("real_wagedata_n=", n))

evalMMDmat<-sapply(1:length(resDRF), function(b)  resDRF[[b]]$evalMMD)
sdMMDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalMMDmat[j,])  )
ylim=c(min(rowMeans(evalMMDmat) - 1.96*sdMMDmat), max(rowMeans(evalMMDmat) + 1.96*sdMMDmat ) ) 
plot(rowMeans(evalMMDmat) - 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen", main="MMD Loss DRF", ylim=ylim, 
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmat) + 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmat), type="l", cex=0.8, col="darkblue")


evalNPLDmat<-sapply(1:length(resDRF), function(b)  resDRF[[b]]$evalNPLD)
sdNPLDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalNPLDmat[j,])  )
ylim=c(min(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat), max(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat ) )
plot(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen", main="NPLD Loss DRF", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalNPLDmat), type="l", cex=0.8, col="darkblue")


evalMMDmat<-sapply(1:length(resDRF_native), function(b)  resDRF_native[[b]]$evalMMD)
sdMMDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalMMDmat[j,])  )
ylim=c(min(rowMeans(evalMMDmat) - 1.96*sdMMDmat), max(rowMeans(evalMMDmat) + 1.96*sdMMDmat ) ) 
plot(rowMeans(evalMMDmat) - 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen", main="MMD Loss DRF native", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmat) + 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmat), type="l", cex=0.8, col="darkblue")


evalNPLDmat<-sapply(1:length(resDRF_native), function(b)  resDRF_native[[b]]$evalNPLD)
sdNPLDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalNPLDmat[j,])  )
ylim=c(min(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat), max(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat ) )
plot(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen", main="NPLD Loss DRF native", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalNPLDmat), type="l", cex=0.8, col="darkblue")


evalMMDmat<-sapply(1:length(ressobolMDA), function(b)  ressobolMDA[[b]]$evalMMD)
sdMMDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalMMDmat[j,])  )
ylim=c(min(rowMeans(evalMMDmat) - 1.96*sdMMDmat), max(rowMeans(evalMMDmat) + 1.96*sdMMDmat ) ) 
plot(rowMeans(evalMMDmat) - 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen", main="MMD Loss Sobol MDA", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmat) + 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmat), type="l", cex=0.8, col="darkblue")


evalNPLDmat<-sapply(1:length(ressobolMDA), function(b)  ressobolMDA[[b]]$evalNPLD)
sdNPLDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalNPLDmat[j,])  )
ylim=c(min(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat), max(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat ) )
plot(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen", main="NPLD Loss Sobol MDA", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalNPLDmat), type="l", cex=0.8, col="darkblue")


##########################################
### Start Evaluation of birth data #######
##########################################




set.seed(10)
n<-2000
ntest <- round(1/2*n)
num.trees<-500

set.seed(10)


tmp<-genData(dataset = "real_birthdata2", n = n)

#tmp<-genData(dataset = "motivatingexample", n = n)


X<-tmp$X
Y<-as.matrix(tmp$y)

if (is.null(colnames(X))){
  colnames(X) <- paste0("X",1:ncol(X))
}


# remove redundant variables
X$race_mother_white <- NULL
X$race_father_white <- NULL
X$marital_status_mother_unmarried <- NULL
X$gender_M <- NULL
X$`delivery_method_C-section` <- NULL


Xtest <- X[(round(n - ntest) + 1):n, , drop = F]
Ytest <- Y[(round(n - ntest) + 1):n, , drop = F]
#
X <- X[1:round(n - ntest), , drop = F]
Y <- Y[1:round(n - ntest), , drop = F]



for (b in 1:B){


VIMPDRF<-compute_drf_vimp(X, Y, X_test = NULL, num.trees = 500, silent = FALSE)
(VIMPDRF<-sort(VIMPDRF))

resDRF[[b]]<-evalfinal(VIMPDRF, X ,Y ,Xtest, Ytest, metrics=c("MMD"), num.trees=500 )



# vimp drf
DRF <- drf(X, Y, num.trees = 500)
vimp_drf <- variable_importance(DRF)
vimp_drf<-sort(c(vimp_drf))
names(vimp_drf)<- colnames(X)

resDRF_native[[b]]<-evalfinal(vimp_drf, X ,Y ,Xtest, Ytest, metrics=c("MMD"), num.trees=500 )




# Sobol MDA
sobolMDAj <- list()
for (j in 1:ncol(Y)){
  XY <- as.data.frame(cbind(X, Y[, j]))
  colnames(XY) <- c(paste('X', 1:(ncol(XY)-1), sep=''), 'Y')
  num.trees <- 500
  clock <- Sys.time()
  forest <- sobolMDA::ranger(Y ~., data = XY, num.trees = num.trees, importance = 'sobolMDA')
  duration <- Sys.time() - clock
  print(duration)
  forest$r.squared
  sobolMDAj[[j]] <- forest$variable.importance
}

sobolMDA<-colMeans(do.call(rbind,sobolMDAj))


ressobolMDA[[b]]<-evalfinal(sobolMDA, X ,Y ,Xtest, Ytest, metrics=c("MMD", "NPLD"), num.trees=500 )
}

save(resDRF, resDRF_native, ressobolMDA, X, Y, Xtest, Ytest, n, ntest, file=paste0("real_birthdata2_n=", n))


evalMMDmat<-sapply(1:length(resDRF), function(b)  resDRF[[b]]$evalMMD)
sdMMDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalMMDmat[j,])  )
ylim=c(min(rowMeans(evalMMDmat) - 1.96*sdMMDmat), max(rowMeans(evalMMDmat) + 1.96*sdMMDmat ) ) 
plot(rowMeans(evalMMDmat) - 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen", main="MMD Loss DRF", ylim=ylim, 
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmat) + 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmat), type="l", cex=0.8, col="darkblue")


evalNPLDmat<-sapply(1:length(resDRF), function(b)  resDRF[[b]]$evalNPLD)
sdNPLDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalNPLDmat[j,])  )
ylim=c(min(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat), max(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat ) )
plot(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen", main="NPLD Loss DRF", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalNPLDmat), type="l", cex=0.8, col="darkblue")


evalMMDmat<-sapply(1:length(resDRF_native), function(b)  resDRF_native[[b]]$evalMMD)
sdMMDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalMMDmat[j,])  )
ylim=c(min(rowMeans(evalMMDmat) - 1.96*sdMMDmat), max(rowMeans(evalMMDmat) + 1.96*sdMMDmat ) ) 
plot(rowMeans(evalMMDmat) - 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen", main="MMD Loss DRF native", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmat) + 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmat), type="l", cex=0.8, col="darkblue")


evalNPLDmat<-sapply(1:length(resDRF_native), function(b)  resDRF_native[[b]]$evalNPLD)
sdNPLDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalNPLDmat[j,])  )
ylim=c(min(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat), max(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat ) )
plot(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen", main="NPLD Loss DRF native", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalNPLDmat), type="l", cex=0.8, col="darkblue")


evalMMDmat<-sapply(1:length(ressobolMDA), function(b)  ressobolMDA[[b]]$evalMMD)
sdMMDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalMMDmat[j,])  )
ylim=c(min(rowMeans(evalMMDmat) - 1.96*sdMMDmat), max(rowMeans(evalMMDmat) + 1.96*sdMMDmat ) ) 
plot(rowMeans(evalMMDmat) - 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen", main="MMD Loss Sobol MDA", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmat) + 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmat), type="l", cex=0.8, col="darkblue")


evalNPLDmat<-sapply(1:length(ressobolMDA), function(b)  ressobolMDA[[b]]$evalNPLD)
sdNPLDmat<-sapply(1:nrow(evalMMDmat), function(j) sd(evalNPLDmat[j,])  )
ylim=c(min(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat), max(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat ) )
plot(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen", main="NPLD Loss Sobol MDA", ylim=ylim,
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalNPLDmat) + 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen")
lines(rowMeans(evalNPLDmat), type="l", cex=0.8, col="darkblue")












