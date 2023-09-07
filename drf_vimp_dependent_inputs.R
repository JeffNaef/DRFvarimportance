

# load packages
library(Matrix)
library(kernlab)
library(drf)
library(MASS)

source('compute_drf_vimp.R')

# define data parameters
p <- 50 # input dimension
n_sample <- 1000 # sample size
rho <- 0.5 # correlation coefficient
beta <- c(1, 1, rep(0, p-2)) # coefficients of linear mean function (X1 and X2 have non-null coefs)
alpha <- c(0, 0, rep(2, 3), rep(0, p-5)) # coefficients of conditional variance (X3, X4, and X5 have non-null coefs)
cov_matrix <- matrix(rho, nrow = p, ncol = p) # covariance matrix
diag(cov_matrix) <- rep(1, p)

# compute drf vimp
niter <- 10 # number of Monte-Carlo repetitions
vimpMC <- sapply(1:niter, function(iter){
  print(paste0('Running iteration number: ', iter, '...'))
  # sample training and testing datasets from Gaussian distributions
  X <- mvrnorm(n = n_sample, mu = rep(0, p), Sigma = cov_matrix)
  X_test <- mvrnorm(n = n_sample, mu = rep(0, p), Sigma = cov_matrix)
  Y <- rnorm(n_sample, mean = X%*%beta, sd = abs(X)%*%alpha)
  # compute vimp
  vimp <- compute_drf_vimp(X, Y, X_test = NULL, num.trees = 500, silent = T)
  return(vimp)
})
vimp <- apply(vimpMC, 1, mean)
vimp_sd <- apply(vimpMC, 1, sd)

# print importance results and computation time
df_vimp <- data.frame(cbind(paste0('X', 1:p), round(vimp, 4),  round(vimp_sd, 4)))
colnames(df_vimp) <- c('Variable', 'vimp', 'std')
df_vimp <- df_vimp[order(-vimp),]
print(df_vimp)

# # check output sampling with Monte-Carlo estimates
# Y_test <- rnorm(n_sample, mean = X_test%*%beta, sd = abs(X_test)%*%alpha)
# Y_check <- sapply(1:10000, function(k){rnorm(n_sample, mean = X%*%beta, sd = abs(X)%*%alpha)})
# print(head(cbind(X[,1]+X[,2], apply(Y_check, 1, mean)), 10)) # conditional mean
# print(head(cbind(2*(abs(X[,3])+abs(X[,4])+abs(X[,5])), apply(Y_check, 1, sd)), 10)) # conditional variance

