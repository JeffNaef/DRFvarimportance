
library(kernlab)
library(drf)
library(Matrix)
# Load necessary libraries
library(mvtnorm) # for generating multivariate normal random variables
library(ks)


library(dplyr)
library(kableExtra)
library(copula)

source("drfnew_v2.R")
#source("applications")
source("genData.R")
source("evaluation.R")

set.seed(10)

### Continue with more fancy examples!!!

n<-500
ntest<-round(n*0.1)
num.trees<-2000

##### Continue with the collection of the DRF datasets ######
## Discuss with Julie: Could be super interesting to use a medical dataset with
## patients, whereby X is the patients characteristics and Y is a stream of measurements
##(though would probably need Gaussian processes in a first step)



## Step 2: Do Analysis

### 2.a) if the dataset is synthetic, we can check the correct variable ordering
evalsynthetic(dataset="motivatingexample", L=10, n=n, B=1, p=2, num.trees = num.trees)



## Step 3: Check quantile error when using only X_1 instead of X_1 and X_2.


tmp<-genData(dataset = "motivatingexample", n = 2*n, p = 2)


Xtrain<-tmp$X[1:n,]
Xtest<-tmp$X[(n+1):(2*n),]
Ytrain<-as.matrix(tmp$y[1:n])
Ytest <-as.matrix(tmp$y[(n+1):(2*n)])
colnames(Xtrain) <- paste0("X",1:ncol(Xtrain))


# target quantiles: 0.025, 0.5, 0.975

DRFall<-drf(X=Xtrain, Y=Ytrain, num.trees=num.trees)
quantpredictall<-predict(DRFall, newdata=Xtest, functional="quantile",quantiles=c(0.025, 0.5, 0.975))

# only with X_1
DRFX1<-drf(X=Xtrain[,1,drop=F], Y=Ytrain, num.trees=num.trees)
quantpredictX1<-predict(DRFX1, newdata=Xtest[,1, drop=F], functional="quantile",quantiles=c(0.025, 0.5, 0.975))


### Compare to actual quantiles:
i<-0
Lossall<-matrix(NaN, nrow=3)
LossX1<-matrix(NaN, nrow=3)

for (tau in c(0.025, 0.5, 0.975)){
  
  i<-i+1
  
  Losslist<-lapply(1:nrow(Xtest),  function(j)  {
    
    truth=qnorm(tau, mean= 0.8*(Xtest[j,1] > 0), sd = sqrt( 1 + (Xtest[j,2] > 0) ) )
    
    eall<-quantpredictall$quantile[,,i][j]
    eX1 <-quantpredictX1$quantile[,,i][j]
    
    
    return( list(Lossall= (truth - eall)*tau*(truth >= eall)+
                   (eall- truth )*(1-tau)*(truth < eall),
              LossX1 = (truth - eX1)*tau*(truth >= eX1)+
                (eX1- truth )*(1-tau)*(truth < eX1)))
    
    
  }  )
  
  Lossall[i]<-mean(sapply(1:nrow(Xtest), function(j)  Losslist[[j]]$Lossall ))
  LossX1[i]<-mean(sapply(1:nrow(Xtest), function(j)  Losslist[[j]]$LossX1 ))
}

Lossall
LossX1

(LossX1-Lossall)/Lossall

