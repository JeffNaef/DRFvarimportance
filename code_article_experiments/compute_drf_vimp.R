
#' Variable importance for Distributional Random Forests
#'
#' @param X Matrix with input training data.
#' @param Y Matrix with output training data.
#' @param X_test Matrix with input testing data. If NULL, out-of-bag estimates are used.
#' @param num.trees Number of trees to fit DRF. Default value is 500 trees.
#' @param silent If FALSE, print variable iteration number, otherwise nothing is print. Default is FALSE.
#'
#' @return The list of importance values for all input variables.
#' @export
#'
#' @examples
compute_drf_vimp <- function(X, Y, X_test = NULL, num.trees = 500, silent = FALSE){
  
  # fit initial DRF
  bandwidth_Y <- drf:::medianHeuristic(Y)
  k_Y <- rbfdot(sigma = bandwidth_Y)
  K <- matrix(kernelMatrix(k_Y, Y, Y), ncol = nrow(Y))
  DRF <- drf(X, Y, num.trees = num.trees)
  wall <- predict(DRF, X_test)$weights
  wall <- matrix(wall, ncol = ncol(wall))
  
  # compute normalization constant
  wbar <- colMeans(wall)
  wall_wbar <- sweep(wall, 2, wbar, "-")
  I0 <- as.numeric(sum(diag(wall_wbar %*% K %*% t(wall_wbar))))
  
  # compute drf importance dropping variables one by one
  I <- sapply(1:ncol(X), function(j) {
    if (!silent){print(paste0('Running importance for variable X', j, '...'))}
    DRFj <- drf(X = X[, -j, drop=F], Y = Y, num.trees = num.trees)
    wj <- predict(DRFj, X_test[, -j])$weights
    wj <- matrix(wj, ncol = ncol(wj))
    Ij <- sum(diag((wj - wall) %*% K %*% t(wj - wall)))/I0
    return(Ij)
  })
  
  # compute retraining bias
  DRF0 <- drf(X = X, Y = Y, num.trees = num.trees)
  w0 = predict(DRF0, X_test)$weights
  w0 <- matrix(w0, ncol = ncol(w0))
  vimp0 <- sum(diag((w0 - wall) %*% K %*% t(w0 - wall)))/I0
  
  # compute final importance (remove bias & truncate negative values)
  vimp <- sapply(I - vimp0, function(x){max(0,x)})
  
  names(vimp)<-colnames(X)
  
  return(vimp)
  
}
