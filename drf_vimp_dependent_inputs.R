

# load packages
library(Matrix)
library(kernlab)
library(drf)
library(MASS)

# define data parameters
p <- 50 # input dimension
n_sample <- 1000 # sample size
rho <- 0.5 # correlation coefficient
beta <- c(1, 1, rep(0, p-2)) # coefficients of linear mean function (X1 and X2 have non-null coefs)
alpha <- c(0, 0, rep(2, 3), rep(0, p-5)) # coefficients of conditional variance (X3, X4, and X5 have non-null coefs)
cov_matrix <- matrix(rho, nrow = p, ncol = p) # covariance matrix
diag(cov_matrix) <- rep(1, p)

# sample training and testing datasets from Gaussian distributions
X <- mvrnorm(n = n_sample, mu = rep(0, p), Sigma = cov_matrix)
X_test <- mvrnorm(n = n_sample, mu = rep(0, p), Sigma = cov_matrix)
Y <- rnorm(n_sample, mean = X%*%beta, sd = abs(X)%*%alpha)
Y_test <- rnorm(n_sample, mean = X_test%*%beta, sd = abs(X_test)%*%alpha)

# # check output sampling with Monte-Carlo estimates
# Y_check <- sapply(1:10000, function(k){rnorm(n_sample, mean = X%*%beta, sd = abs(X)%*%alpha)})
# print(head(cbind(X[,1]+X[,2], apply(Y_check, 1, mean)), 10)) # conditional mean
# print(head(cbind(2*(abs(X[,3])+abs(X[,4])+abs(X[,5])), apply(Y_check, 1, sd)), 10)) # conditional variance

# fit initial DRF
num_trees <- 500
bandwidth_Y <- drf:::medianHeuristic(Y_test)
k_Y <- rbfdot(sigma = bandwidth_Y)
K <- kernelMatrix(k_Y, Y, y = Y)
DRF <- drf(X, Y, num.trees = num_trees)
wall <- predict(DRF, x = X_test)$weights

# compute normalization constant
wbar <- colMeans(wall)
wall_wbar <- sweep(wall, 2, wbar, "-")
I0 <- as.numeric(sum(diag(wall_wbar %*% K %*% t(wall_wbar))))

# compute importance for all variables
time_initial <- Sys.time()
I <- sapply(1:p, function(j) {
  print(paste0('Running importance for variable X', j, '...'))
  DRFj <- drf(X = X[, -j, drop=F], Y = Y, num.trees = num_trees) # refit drf dropping variables one by one
  DRFpredj <- predict(DRFj, x = X_test[, -j])
  wj <- DRFpredj$weights
  Ij <- sum(diag((wj - wall) %*% K %*% t(wj - wall)))/I0
  return(Ij)
})

# compute retraining bias
DRF0 <- drf(X = X, Y = Y, num.trees = num_trees)
DRFpred0 = predict(DRF0, x = X_test)
w0 <- DRFpred0$weights
vimp0 <- sum(diag((w0 - wall) %*% K %*% t(w0 - wall)))/I0

# compute final importance (remove bias & truncate negative values)
vimp <- sapply(I - vimp0, function(x){max(0,x)})
time_final <- Sys.time()

# print importance results and computation time
df_vimp <- data.frame(cbind(paste0('X', 1:p), round(vimp, 3)))
colnames(df_vimp) <- c('Variable', 'vimp')
df_vimp <- df_vimp[order(-vimp),]
print(df_vimp)
print(paste0('Computation time (min): ', time_final - time_initial))
