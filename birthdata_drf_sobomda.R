
# load packages
library(Matrix)
library(kernlab)
library(drf)
library(MASS)
library(sobolMDA)

source('compute_drf_vimp.R')

# sample data
load("applications/births_data/data/datasets/births_benchmark.Rdata")
X <- as.data.frame(X)
dim(Y)
dim(X)
# remove redundant variables
X$race_mother_white <- NULL
X$race_father_white <- NULL
X$marital_status_mother_unmarried <- NULL
X$gender_M <- NULL
X$`delivery_method_C-section` <- NULL

nsample <- 2000
p <- ncol(X)
ind <- 1:nrow(X)
ind_train <- sample(ind, nsample, replace = F)
X_train <- X[ind_train,]
Y_train <- Y[ind_train, ]
ind_test <- sample(setdiff(ind, ind_train), nsample, replace = F)
X_test <- X[ind_test,]
Y_test <- Y[ind_test, ]


# compute drf vimp
vimp <- compute_drf_vimp(X_train, Y_train)

# vimp drf
bandwidth_Y <- drf:::medianHeuristic(Y_train)
k_Y <- rbfdot(sigma = bandwidth_Y)
K <- kernelMatrix(k_Y, Y_train, Y_train)
DRF <- drf(X_train, Y_train, num.trees = 500)
vimp_drf <- variable_importance(DRF)

# compute SobolMDA
sobolMDA <- list()
for (j in 1:ncol(Y)){
  XY <- as.data.frame(cbind(X_train, Y_train[, j]))
  colnames(XY) <- c(paste('X', 1:(ncol(XY)-1), sep=''), 'Y')
  num.trees <- 300
  clock <- Sys.time()
  forest <- sobolMDA::ranger(Y ~., data = XY, num.trees = num.trees, importance = 'sobolMDA')
  duration <- Sys.time() - clock
  print(duration)
  forest$r.squared
  sobolMDA[[j]] <- forest$variable.importance
}

# print results
df_vimp <- data.frame(colnames(X), round(vimp, 4), round(vimp_drf, 4), round(sobolMDA[[1]], 4), round(sobolMDA[[2]], 4))
colnames(df_vimp) <- c('Variable', 'vimp', 'vimp-drf', 'sobolMDA-Y1', 'sobolMDA-Y2')
df_vimp <- df_vimp[order(-vimp),]
print(df_vimp)


