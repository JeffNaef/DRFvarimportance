

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





## Step 2: Do Analysis

n<-500
ntest<-round(n*0.1)
num.trees<-2000

set.seed(10)
### 2.a) if the dataset is synthetic, we can check the correct variable ordering
evalsynthetic(dataset="copulasynthetic", L=10, n=n, B=1, p=10, num.trees = num.trees, MRF=T)

