

# load packages
library(Matrix)
library(kernlab)
library(drf)
library(MASS)
library(sobolMDA) # package can be found here https://gitlab.com/drti/sobolmda (in folder release)
library(parallel)
library(xtable)

source('compute_drf_vimp.R')


# define data parameters
p <- 10 # input dimension
nsample <- 2000 # sample size
rho <- 0.5 # correlation coefficient
beta <- c(2, 1, rep(0, p-2)) # coefficients of linear mean function (X1 and X2 have non-null coefs)
alpha <- c(0, 0, rep(2, 3), rep(0, p-5)) # coefficients of conditional variance (X3, X4, and X5 have non-null coefs)
cov_matrix <- matrix(rho, nrow = p, ncol = p) # covariance matrix
diag(cov_matrix) <- rep(1, p)
cov_matrix[1, p] <- cov_matrix[p, 1] <- 0.9


# define algorithm parameters
num.trees <- 500
niter <- 10


# define cluster
ncores <- 6
clust <- makeCluster(ncores)
clusterExport(clust, list("p", "nsample", "beta", "alpha", "cov_matrix", "num.trees",
                          "drf", "compute_drf_vimp", "variable_importance", "ranger",
                          "mvrnorm", "rbfdot", "diag", "kernelMatrix"))


# run vimp
time_start <- Sys.time()
drf_vimp_MC <- parLapply(clust, 1:niter, function(iter){

  # sample data
  X <- mvrnorm(n = nsample, mu = rep(0, p), Sigma = cov_matrix)
  Y <- matrix(rnorm(nsample, mean = X%*%beta, sd = abs(X)%*%alpha), ncol = 1)
  
  # DRF vimp
  vimp <- compute_drf_vimp(X, Y, silent = F)
  
  # vimp drf
  DRF <- drf(X, Y, num.trees = num.trees)
  vimp_drf <- variable_importance(DRF)
  
  # sobolMDA
  XY <- as.data.frame(cbind(X, Y))
  colnames(XY) <- c(paste('X', 1:(ncol(XY)-1), sep=''), 'Y')
  forest <- sobolMDA::ranger(Y ~., data = XY, num.trees = num.trees, importance = 'sobolMDA')
  vimp_sobol <- forest$variable.importance
  forest <- sobolMDA::ranger(Y ~., data = XY, num.trees = num.trees, importance = 'permutation')
  vimp_MDA <- forest$variable.importance
  
  return(cbind(vimp, vimp_drf, vimp_MDA, vimp_sobol))
  
})
duration <- Sys.time() - time_start
print(duration)

vimp_In <- sapply(drf_vimp_MC, function(df){
  df[, 1]
})
vimp_drf <- sapply(drf_vimp_MC, function(df){
  df[, 2]
})
vimp_MDA <- sapply(drf_vimp_MC, function(df){
  df[, 3]
})
vimp_sobol <- sapply(drf_vimp_MC, function(df){
  df[, 4]
})

# save results
save.image(file = paste0('drf_vimp_xp1.Rdata'))

# print results
cbind(apply(vimp_In, 1, mean), apply(vimp_drf, 1, mean), apply(vimp_MDA, 1, mean), apply(vimp_sobol, 1, mean))
cbind(apply(vimp_In, 1, sd), apply(vimp_drf, 1, sd), apply(vimp_MDA, 1, sd), apply(vimp_sobol, 1, sd))/sqrt(niter)

df_In <- cbind(apply(vimp_In, 1, mean), apply(vimp_In, 1, sd)/sqrt(niter))
df_In <- df_In[order(-df_In[,1]),]
round(df_In, 4)
xtable(df_In, digits = 4)

df_drf <- cbind(apply(vimp_drf, 1, mean), apply(vimp_drf, 1, sd)/sqrt(niter))
df_drf <- df_drf[order(-df_drf[,1]),]
round(df_drf, 4)
xtable(df_drf, digits = 4)

df_MDA <- cbind(apply(vimp_MDA, 1, mean), apply(vimp_MDA, 1, sd)/sqrt(niter))
df_MDA <- df_MDA[order(-df_MDA[,1]),]
round(df_MDA, 2)
xtable(df_MDA, digits = 4)

df_sobol <- cbind(apply(vimp_sobol, 1, mean), apply(vimp_sobol, 1, sd)/sqrt(niter))
df_sobol <- df_sobol[order(-df_sobol[,1]),,drop=F]
round(df_sobol, 4)
xtable(df_sobol, digits = 4)
