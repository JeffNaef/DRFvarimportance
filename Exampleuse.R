
source("drfnew_v2.R")
library(kernlab)
library(drf)

n<-500
p<-4




## Example 1
## In this example X1 generates a covariate shift in Y
x <- runif(n,-1,1)
Y <- as.matrix(rnorm(n, mean = 0.8*(x > 0)))
X <- matrix(runif(n * (p-1)), ncol = p-1)
X <- cbind(x,X)
colnames(X)<-paste0("X",1:p)


res<-drfwithVI(X, Y, B=1, num.trees=3000, num.features=10, min.node.size=15)# sample.splitting = T, ntest=100

# Variable importance 
res$VI



## Example 2
## In this example X1 generates a quantile shift in Y (from light tailed to heavy- tailed)
x <- runif(n,-1,1)
Y <- ifelse(x >= 0, rt(n = n, 3), rnorm(n, 1, 1))
X <- matrix(runif(n * (p-1)), ncol = p-1)
X <- cbind(x,X)
colnames(X)<-paste0("X",1:p)


res<-drfwithVI(X, Y, B=1, num.trees=3000, num.features=10, min.node.size=15)# sample.splitting = T, ntest=100

# Variable importance 
res$VI

