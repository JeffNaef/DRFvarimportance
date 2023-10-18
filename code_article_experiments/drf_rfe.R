

# load packages
library(Matrix)
library(kernlab)
library(MASS)
library(foreign)
library(parallel)
library(drf)


# define utility functions
source('compute_drf_vimp.R')
source('load_mtr_datasets.R')
compute_pinball_loss <- function(Y_test, Y_quantiles, level){
  x <- Y_test - Y_quantiles
  level*mean(abs(x[x > 0])) + (1-level)*mean(abs(x[x < 0]))
}
compute_mmd_loss <- function(Y_train, Y_test, weights){
  # Y_train <- scale(Y_train)
  # Y_test <- scale(Y_test)
  bandwidth_Y <- (1/drf:::medianHeuristic(Y_train))^2
  k_Y <- rbfdot(sigma = bandwidth_Y)
  K_train <- matrix(kernelMatrix(k_Y, Y_train, Y_train), ncol = nrow(Y_train))
  K_cross <- matrix(kernelMatrix(k_Y, Y_test, Y_train), ncol = nrow(Y_train))
  weights <- matrix(weights, ncol = ncol(weights))
  t1 <- diag(weights%*%K_train%*%t(weights))
  t2 <- diag(K_cross%*%t(weights))
  mmd_loss <- mean(t1) - 2*mean(t2)
  mmd_loss
}
compute_rfe <- function(iter, algo){
  
  print(paste0("Iteration number: ", iter))
  
  X <- data$X
  Y <- data$Y
  p <- data$p
  
  # subsampling for large datasets
  nsample <- nrow(X)
  nmax <- 2000
  if (nsample > nmax){
    X <- X[subsampling_index[[iter]], ]
    Y <- Y[subsampling_index[[iter]], ]
  }

  # split data in two parts
  nsample <- round(nrow(X)/2)
  ind_train <- sample(1:nrow(X), nsample, replace = F)
  X_train <- X[ind_train,]
  Y_train <- Y[ind_train,]
  ind_test <- setdiff(1:nrow(X), ind_train)
  X_test <- X[ind_test,]
  Y_test <- Y[ind_test,]

  # run RFE
  variable_list <- c()
  mmd_list <- c()
  X_train_temp <- X_train
  X_test_temp <- X_test
  for (j in 1:p){
    DRFj <- drf(X_train_temp, Y_train, num.trees = 500)
    weights <- predict(DRFj, X_test_temp)$weights
    mmd <- compute_mmd_loss(Y_train, Y_test, weights)
    mmd_list <- c(mmd, mmd_list)
    if (algo == "drf_freq"){
      vimpj <- variable_importance(DRFj)
    }
    if (algo == "vimp"){
      vimpj <- compute_drf_vimp(X_train_temp, Y_train, silent = T)
    }
    ind_min <- which.min(vimpj)
    # variable_list <- c(colnames(X_train_temp)[ind_min], variable_list)
    X_train_temp <- X_train_temp[, -ind_min, drop=F]
    X_test_temp <- X_test_temp[, -ind_min, drop=F]
  }

  return(mmd_list)
  
}


# load dataset
args <- commandArgs(trailingOnly=TRUE)
dataset <- args[1]
niter <- args[2]
scaling <- args[3]
print(c(dataset, niter, scaling))
data <- load_mtr_datasets(dataset)
if (scaling){
  print('Scaling output')
  data$Y <- scale(data$Y)
}


# subsampling for large datasets
nsample <- nrow(data$X)
nmax <- 2000
subsampling_index <- list()
if (nsample > nmax){
  for (iter in 1:niter){
    subsampling_index <- lapply(1:niter, function(iter){
      sample(1:nsample, nmax, replace = F)
    })
  }
}


# run RFE
ncores <- detectCores()
clust <- makeCluster(ncores)
clusterExport(clust, list("data", "dataset", "compute_rfe", "predict", "subsampling_index",
                          "drf", "rbfdot", "compute_mmd_loss", "compute_drf_vimp", "scaling",
                          "variable_importance", "kernelMatrix", "diag"))
# basic DRF importance
mmd_rfe_drf <- parLapply(clust, 1:niter, compute_rfe, algo = "drf_freq")
mmd_rfe_drf <- do.call('cbind', mmd_rfe_drf)
# vimp
mmd_rfe_vimp <- parLapply(clust, 1:niter, compute_rfe, algo = 'vimp')
mmd_rfe_vimp <- do.call('cbind', mmd_rfe_vimp)


# save results
save.image(file = paste0('mtr_results/rfe_results_', dataset, '.Rdata'))


