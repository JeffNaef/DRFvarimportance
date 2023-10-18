
# load packages
library(Matrix)
library(kernlab)
library(drf)
library(MASS)
library(sobolMDA) # package can be found here https://gitlab.com/drti/sobolmda (in folder release)
library(parallel)
library(xtable)
library(mvtnorm)

source('compute_drf_vimp.R')

# data parameters
xp_num <- 2 # 2 for the bivariate case and 3 for the functional case.
p <- 10
nsample <- 2000
sample_data <- function(xp_num){
  X <- matrix(runif(nsample * p, 0, 1), ncol = p) 
  Y <- NULL
  if (xp_num == 2){
    Y1 <- runif(nsample, min = X[, 1], max = X[, 1] + 1)
    Y2 <- runif(nsample, min = rep(0, nsample), max = X[, 2])
    Y <- cbind(Y1, Y2)
  }
  if (xp_num == 3){
    d <- 30
    t <- seq(-5, 5, length.out = d)
    Y <- t(sapply(1:nsample, function(i){
      lengthscale <- X[i, 2]^2  # lengthscale parameter for the RBF kernel depending on X_2
      K <- kernelMatrix(rbfdot(sigma =lengthscale), t, y = t)
      rmvnorm(1, mean = rep(X[i, 1], d), sigma = K)
    }))
  }
  return(list(X = X, Y = Y))
}

# define algorithm parameters
num.trees <- 500
niter <- 10


# define cluster
ncores <- 6
clust <- makeCluster(ncores)
clusterExport(clust, list("p", "nsample", "sample_data", "num.trees", "xp_num", "rmvnorm",
                          "drf", "compute_drf_vimp", "variable_importance", "ranger",
                          "mvrnorm", "rbfdot", "diag", "kernelMatrix"))


# run vimp
time_start <- Sys.time()
drf_vimp_MC <- parLapply(clust, 1:niter, function(iter){

  # sample data
  data <- sample_data(xp_num)
  X <- data$X
  Y <- data$Y
  
  # DRF vimp
  vimp <- compute_drf_vimp(X, Y, silent = F)
  
  # vimp drf
  DRF <- drf(X, Y, num.trees = num.trees)
  vimp_drf <- variable_importance(DRF)

  return(cbind(vimp, vimp_drf))
  
})

# print results
vimp_In <- sapply(drf_vimp_MC, function(df){
  df[, 1]
})
vimp_drf <- sapply(drf_vimp_MC, function(df){
  df[, 2]
})

df_In <- cbind(apply(vimp_In, 1, mean), apply(vimp_drf, 1, mean), apply(vimp_In, 1, sd)/sqrt(niter), apply(vimp_drf, 1, sd)/sqrt(niter))
rownames(df_In) <- paste0('X', 1:p)
df_In <- df_In[order(-df_In[,1]),]
round(df_In, 4)
xtable(t(df_In[, 1:2]), digits = 2)

