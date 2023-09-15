

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
library(latex2exp)


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
names(sobolMDA)<- colnames(X)


ressobolMDA[[b]]<-evalfinal(sobolMDA, X ,Y ,Xtest, Ytest, metrics=c("MMD","NPLD","MAD"), num.trees=500 )
}

save(resDRF, resDRF_native, ressobolMDA, X,Y, Xtest, Ytest, n, ntest, file=paste0("real_wagedata_n=", n))


####### MMD #####################

evalMMDmatDRF<-sapply(1:length(resDRF), function(b)  resDRF[[b]]$evalMMD)
sdMMDmatDRF<-sapply(1:nrow(evalMMDmatDRF), function(j) sd(evalMMDmatDRF[j,])  )

evalMMDmatDRFnative<-sapply(1:length(resDRF_native), function(b)  resDRF_native[[b]]$evalMMD)
sdMMDmatDRFnative<-sapply(1:nrow(evalMMDmatDRFnative), function(j) sd(evalMMDmatDRFnative[j,])  )


evalMMDmatSobol<-sapply(1:length(ressobolMDA), function(b)  ressobolMDA[[b]]$evalMMD)
sdMMDmatSobol<-sapply(1:nrow(evalMMDmatSobol), function(j) sd(evalMMDmatSobol[j,])  )


ylim1=c(min(rowMeans(evalMMDmatDRF) - 1.96*sdMMDmatDRF), max(rowMeans(evalMMDmatDRF) + 1.96*sdMMDmatDRF ) ) 
ylim2=c(min(rowMeans(evalMMDmatDRFnative) - 1.96*sdMMDmatDRFnative), max(rowMeans(evalMMDmatDRFnative) + 1.96*sdMMDmatDRFnative ) ) 
ylim3=c(min(rowMeans(evalMMDmatSobol) - 1.96*sdMMDmatSobol), max(rowMeans(evalMMDmatSobol) + 1.96*sdMMDmatSobol ) ) 


ylim=c(min( c(ylim1[1] , ylim2[1], ylim3[1]  )), max( c(ylim1[2] , ylim2[2], ylim3[2]  )     ) )


## ours
#TeX(r'(MMD Loss $I_n^{(j)}$)', bold=TRUE)
plot(rowMeans(evalMMDmatDRF) - 1.96*sdMMDmatDRF, type="l", lwd=2, cex=0.8, col="darkgreen", main="MMD loss" , ylim=ylim, 
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmatDRF) + 1.96*sdMMDmatDRF, type="l", lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmatDRF), type="l",lwd=2, cex=0.8, col="darkblue")


#DRF native
#ylim=c(min(rowMeans(evalMMDmatDRFnative) - 1.96*sdMMDmatDRFnative), max(rowMeans(evalMMDmatDRFnative) + 1.96*sdMMDmatDRFnative ) ) 
lines(rowMeans(evalMMDmatDRFnative) - 1.96*sdMMDmatDRFnative, type="l", lty=2,lwd=2, cex=0.8, col="darkgreen")
#plot(rowMeans(evalMMDmat) - 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen", main="MMD Loss DRF native", ylim=ylim,
#     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmatDRFnative) + 1.96*sdMMDmatDRFnative, type="l", lty=2,lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmatDRFnative), type="l", lty=2,lwd=2, cex=0.8, col="darkblue")



lines(rowMeans(evalMMDmatSobol) - 1.96*sdMMDmatSobol, type="l", lty=3,lwd=2, cex=0.8, col="darkgreen")
#ylim=c(min(rowMeans(evalMMDmatSobol) - 1.96*sdMMDmatSobol), max(rowMeans(evalMMDmatSobol) + 1.96*sdMMDmatSobol ) ) 
#plot(rowMeans(evalMMDmatSobol) - 1.96*sdMMDmatSobol, type="l", cex=0.8, col="darkgreen", main="MMD Loss Sobol-MDA", ylim=ylim,
#     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmatSobol) + 1.96*sdMMDmatSobol, type="l", lty=3,lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmatSobol), type="l", lty=3,lwd=2, cex=0.8, col="darkblue")

####### MMD #####################



####### NPLD #####################

evalNPLDmatDRF<-sapply(1:length(resDRF), function(b)  resDRF[[b]]$evalNPLD)
sdNPLDmatDRF<-sapply(1:nrow(evalNPLDmatDRF), function(j) sd(evalNPLDmatDRF[j,])  )

evalNPLDmatDRFnative<-sapply(1:length(resDRF_native), function(b)  resDRF_native[[b]]$evalNPLD)
sdNPLDmatDRFnative<-sapply(1:nrow(evalNPLDmatDRFnative), function(j) sd(evalNPLDmatDRFnative[j,])  )


evalNPLDmatSobol<-sapply(1:length(ressobolMDA), function(b)  ressobolMDA[[b]]$evalNPLD)
sdNPLDmatSobol<-sapply(1:nrow(evalNPLDmatSobol), function(j) sd(evalNPLDmatSobol[j,])  )


ylim1=c(min(rowMeans(evalNPLDmatDRF) - 1.96*sdNPLDmatDRF), max(rowMeans(evalNPLDmatDRF) + 1.96*sdNPLDmatDRF ) ) 
ylim2=c(min(rowMeans(evalNPLDmatDRFnative) - 1.96*sdNPLDmatDRFnative), max(rowMeans(evalNPLDmatDRFnative) + 1.96*sdNPLDmatDRFnative ) ) 
ylim3=c(min(rowMeans(evalNPLDmatSobol) - 1.96*sdNPLDmatSobol), max(rowMeans(evalNPLDmatSobol) + 1.96*sdNPLDmatSobol ) ) 


ylim=c(min( c(ylim1[1] , ylim2[1], ylim3[1]  )), max( c(ylim1[2] , ylim2[2], ylim3[2]  )     ) )


## ours
#TeX(r'(NPLD Loss $I_n^{(j)}$)', bold=TRUE)
plot(rowMeans(evalNPLDmatDRF) - 1.96*sdNPLDmatDRF, type="l", lwd=2, cex=0.8, col="darkgreen", main="NPLD loss" , ylim=ylim, 
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalNPLDmatDRF) + 1.96*sdNPLDmatDRF, type="l", lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalNPLDmatDRF), type="l",lwd=2, cex=0.8, col="darkblue")


#DRF native
#ylim=c(min(rowMeans(evalNPLDmatDRFnative) - 1.96*sdNPLDmatDRFnative), max(rowMeans(evalNPLDmatDRFnative) + 1.96*sdNPLDmatDRFnative ) ) 
lines(rowMeans(evalNPLDmatDRFnative) - 1.96*sdNPLDmatDRFnative, type="l", lty=2,lwd=2, cex=0.8, col="darkgreen")
#plot(rowMeans(evalNPLDmat) - 1.96*sdNPLDmat, type="l", cex=0.8, col="darkgreen", main="NPLD Loss DRF native", ylim=ylim,
#     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalNPLDmatDRFnative) + 1.96*sdNPLDmatDRFnative, type="l", lty=2,lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalNPLDmatDRFnative), type="l", lty=2,lwd=2, cex=0.8, col="darkblue")



lines(rowMeans(evalNPLDmatSobol) - 1.96*sdNPLDmatSobol, type="l", lty=3,lwd=2, cex=0.8, col="darkgreen")
#ylim=c(min(rowMeans(evalNPLDmatSobol) - 1.96*sdNPLDmatSobol), max(rowMeans(evalNPLDmatSobol) + 1.96*sdNPLDmatSobol ) ) 
#plot(rowMeans(evalNPLDmatSobol) - 1.96*sdNPLDmatSobol, type="l", cex=0.8, col="darkgreen", main="NPLD Loss Sobol-MDA", ylim=ylim,
#     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalNPLDmatSobol) + 1.96*sdNPLDmatSobol, type="l", lty=3,lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalNPLDmatSobol), type="l", lty=3,lwd=2, cex=0.8, col="darkblue")

####### NPLD #####################


####### MAD #####################

evalMADmatDRF<-sapply(1:length(resDRF), function(b)  resDRF[[b]]$evalMAD)
sdMADmatDRF<-sapply(1:nrow(evalMADmatDRF), function(j) sd(evalMADmatDRF[j,])  )

evalMADmatDRFnative<-sapply(1:length(resDRF_native), function(b)  resDRF_native[[b]]$evalMAD)
sdMADmatDRFnative<-sapply(1:nrow(evalMADmatDRFnative), function(j) sd(evalMADmatDRFnative[j,])  )


evalMADmatSobol<-sapply(1:length(ressobolMDA), function(b)  ressobolMDA[[b]]$evalMAD)
sdMADmatSobol<-sapply(1:nrow(evalMADmatSobol), function(j) sd(evalMADmatSobol[j,])  )


ylim1=c(min(rowMeans(evalMADmatDRF) - 1.96*sdMADmatDRF), max(rowMeans(evalMADmatDRF) + 1.96*sdMADmatDRF ) ) 
ylim2=c(min(rowMeans(evalMADmatDRFnative) - 1.96*sdMADmatDRFnative), max(rowMeans(evalMADmatDRFnative) + 1.96*sdMADmatDRFnative ) ) 
ylim3=c(min(rowMeans(evalMADmatSobol) - 1.96*sdMADmatSobol), max(rowMeans(evalMADmatSobol) + 1.96*sdMADmatSobol ) ) 


ylim=c(min( c(ylim1[1] , ylim2[1], ylim3[1]  )), max( c(ylim1[2] , ylim2[2], ylim3[2]  )     ) )


## ours
#TeX(r'(MAD Loss $I_n^{(j)}$)', bold=TRUE)
plot(rowMeans(evalMADmatDRF) - 1.96*sdMADmatDRF, type="l", lwd=2, cex=0.8, col="darkgreen", main="MAD loss" , ylim=ylim, 
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMADmatDRF) + 1.96*sdMADmatDRF, type="l", lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalMADmatDRF), type="l",lwd=2, cex=0.8, col="darkblue")


#DRF native
#ylim=c(min(rowMeans(evalMADmatDRFnative) - 1.96*sdMADmatDRFnative), max(rowMeans(evalMADmatDRFnative) + 1.96*sdMADmatDRFnative ) ) 
lines(rowMeans(evalMADmatDRFnative) - 1.96*sdMADmatDRFnative, type="l", lty=2,lwd=2, cex=0.8, col="darkgreen")
#plot(rowMeans(evalMADmat) - 1.96*sdMADmat, type="l", cex=0.8, col="darkgreen", main="MAD Loss DRF native", ylim=ylim,
#     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMADmatDRFnative) + 1.96*sdMADmatDRFnative, type="l", lty=2,lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalMADmatDRFnative), type="l", lty=2,lwd=2, cex=0.8, col="darkblue")



lines(rowMeans(evalMADmatSobol) - 1.96*sdMADmatSobol, type="l", lty=3,lwd=2, cex=0.8, col="darkgreen")
#ylim=c(min(rowMeans(evalMADmatSobol) - 1.96*sdMADmatSobol), max(rowMeans(evalMADmatSobol) + 1.96*sdMADmatSobol ) ) 
#plot(rowMeans(evalMADmatSobol) - 1.96*sdMADmatSobol, type="l", cex=0.8, col="darkgreen", main="MAD Loss Sobol-MDA", ylim=ylim,
#     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMADmatSobol) + 1.96*sdMADmatSobol, type="l", lty=3,lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalMADmatSobol), type="l", lty=3,lwd=2, cex=0.8, col="darkblue")

####### MAD #####################




##########################################
### Start Evaluation of birth data #######
##########################################




set.seed(10)
n<-2000
ntest <- round(1/2*n)
num.trees<-500


set.seed(10)

B<-10

resDRF<-list()
resDRF_native<-list()
ressobolMDA<-list()


tmp<-genData(dataset = "real_birthdata2", n = n)

#tmp<-genData(dataset = "motivatingexample", n = n)


X<-tmp$X
Y<-as.matrix(tmp$y)

if (is.null(colnames(X))){
  colnames(X) <- paste0("X",1:ncol(X))
}


# remove redundant variables

X<-X[, !(colnames(X) %in% c("race_mother_white", "race_father_white",  "marital_status_mother_unmarried",
                            "gender_M", "delivery_method_C-section" ) )  ]
# X[,"race_mother_white"] <- NULL
# X[,"race_father_white"] <- NULL
# X[,"marital_status_mother_unmarried"] <- NULL
# X[,"gender_M"] <- NULL
# X[,"delivery_method_C-section"] <- NULL


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
names(sobolMDA)<- colnames(X)


ressobolMDA[[b]]<-evalfinal(sobolMDA, X ,Y ,Xtest, Ytest, metrics=c("MMD"), num.trees=500 )
}

save(resDRF, resDRF_native, ressobolMDA, X, Y, Xtest, Ytest, n, ntest, file=paste0("real_birthdata2_n=", n))


####### MMD #####################

evalMMDmatDRF<-sapply(1:length(resDRF), function(b)  resDRF[[b]]$evalMMD)
sdMMDmatDRF<-sapply(1:nrow(evalMMDmatDRF), function(j) sd(evalMMDmatDRF[j,])  )

evalMMDmatDRFnative<-sapply(1:length(resDRF_native), function(b)  resDRF_native[[b]]$evalMMD)
sdMMDmatDRFnative<-sapply(1:nrow(evalMMDmatDRFnative), function(j) sd(evalMMDmatDRFnative[j,])  )


evalMMDmatSobol<-sapply(1:length(ressobolMDA), function(b)  ressobolMDA[[b]]$evalMMD)
sdMMDmatSobol<-sapply(1:nrow(evalMMDmatSobol), function(j) sd(evalMMDmatSobol[j,])  )


ylim1=c(min(rowMeans(evalMMDmatDRF) - 1.96*sdMMDmatDRF), max(rowMeans(evalMMDmatDRF) + 1.96*sdMMDmatDRF ) ) 
ylim2=c(min(rowMeans(evalMMDmatDRFnative) - 1.96*sdMMDmatDRFnative), max(rowMeans(evalMMDmatDRFnative) + 1.96*sdMMDmatDRFnative ) ) 
ylim3=c(min(rowMeans(evalMMDmatSobol) - 1.96*sdMMDmatSobol), max(rowMeans(evalMMDmatSobol) + 1.96*sdMMDmatSobol ) ) 


ylim=c(min( c(ylim1[1] , ylim2[1], ylim3[1]  )), max( c(ylim1[2] , ylim2[2], ylim3[2]  )     ) )


## ours
#TeX(r'(MMD Loss $I_n^{(j)}$)', bold=TRUE)
plot(rowMeans(evalMMDmatDRF) - 1.96*sdMMDmatDRF, type="l", lwd=2, cex=0.8, col="darkgreen", main="MMD loss" , ylim=ylim, 
     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmatDRF) + 1.96*sdMMDmatDRF, type="l", lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmatDRF), type="l",lwd=2, cex=0.8, col="darkblue")


#DRF native
#ylim=c(min(rowMeans(evalMMDmatDRFnative) - 1.96*sdMMDmatDRFnative), max(rowMeans(evalMMDmatDRFnative) + 1.96*sdMMDmatDRFnative ) ) 
lines(rowMeans(evalMMDmatDRFnative) - 1.96*sdMMDmatDRFnative, type="l", lty=2,lwd=2, cex=0.8, col="darkgreen")
#plot(rowMeans(evalMMDmat) - 1.96*sdMMDmat, type="l", cex=0.8, col="darkgreen", main="MMD Loss DRF native", ylim=ylim,
#     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmatDRFnative) + 1.96*sdMMDmatDRFnative, type="l", lty=2,lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmatDRFnative), type="l", lty=2,lwd=2, cex=0.8, col="darkblue")



lines(rowMeans(evalMMDmatSobol) - 1.96*sdMMDmatSobol, type="l", lty=3,lwd=2, cex=0.8, col="darkgreen")
#ylim=c(min(rowMeans(evalMMDmatSobol) - 1.96*sdMMDmatSobol), max(rowMeans(evalMMDmatSobol) + 1.96*sdMMDmatSobol ) ) 
#plot(rowMeans(evalMMDmatSobol) - 1.96*sdMMDmatSobol, type="l", cex=0.8, col="darkgreen", main="MMD Loss Sobol-MDA", ylim=ylim,
#     xlab="Number of Variables removed", ylab="Values")
lines(rowMeans(evalMMDmatSobol) + 1.96*sdMMDmatSobol, type="l", lty=3,lwd=2, cex=0.8, col="darkgreen")
lines(rowMeans(evalMMDmatSobol), type="l", lty=3,lwd=2, cex=0.8, col="darkblue")

####### MMD #####################

