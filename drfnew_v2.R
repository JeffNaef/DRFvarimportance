
drfwithVI <- function(X, Y, B, sampling = "binomial", sample.splitting=F,...){
  
  args<-list(...)
  
  if (sample.splitting == T) {
    
    ntest<-args$ntest
    
    # Sample Splitting
    Xtest <- X[(round(n - ntest) + 1):n, , drop = F]
    Ytest <- Y[(round(n - ntest) + 1):n, , drop = F]
    #
    X <- X[1:round(n - ntest), , drop = F]
    Y <- Y[1:round(n - ntest), , drop = F]
  } else{
    # No sample splitting
    Xtest <- X
    Ytest <- Y
  }
  
  
  
  bandwidth_Y <- drf:::medianHeuristic(Ytest)
  k_Y <- rbfdot(sigma = bandwidth_Y)
  K <- kernelMatrix(k_Y, Y, y = Y)
  
  
  if (B > 1){
    
    DRF <-
      drfCI(
        X = X,
        Y = Y,
        B = B,
        num.trees = num.trees
      )
    
    # Prediction with all X
    DRFpred = predictdrf(DRF, x = Xtest)
    wall <- DRFpred$weights
    ##
    
    
    I0 <- sapply(1:ncol(Xtest), function(j) {
      # iterate over class 1
      
      ## With CI
      DRFj <-
        drfCI(
          X = X[, -j],
          Y = Y,
          B = B,
          num.trees = num.trees
        )
      
      DRFpredj = predictdrf(DRFj, x = Xtest[, -j])
      wj <- DRFpredj$weights
      val <- sum(diag( (wj - wall) %*% K %*% t(wj - wall) ))
      
      
      
      # Get null distribution if B > 1
      nulldist <- sapply(1:B, function(b) {
        # iterate over class 1
        
        wbj <- DRFpredj$weightsb[[b]]
        wb <- DRFpred$weightsb[[b]]
        
        sum(diag((wb - wall - (wbj - wj)) %*% K %*% t(wb - wall - (wbj - wj))))
      })
      ##
      right_quantile <- quantile(nulldist, 1 - alpha)
      
      
      max(val - unname(right_quantile), 0)
      
      
      
    })
    
  }else{
    
    DRF2<-drf(X,Y, num.trees=num.trees)
    wall<-predict(DRF2,x = Xtest)$weights
    
    I0 <- sapply(1:ncol(Xtest), function(j) {
      # iterate over class 1
      
      ## With CI
      DRFj <-
        drf(
          X = X[, -j],
          Y = Y,
          num.trees = num.trees
        )
      
      DRFpredj = predict(DRFj, x = Xtest[, -j])
      wj <- DRFpredj$weights
      val <- sum(diag( (wj - wall) %*% K %*% t(wj - wall) ))
      
    })
    
    
  }
  
  ### Calculate average weight
  wbar <- colMeans(wall)
  # alternative:
  #wbar<- rep(1/ntest,ntest)
  wall_wbar<-sweep(wall, 2, wbar, "-")
  #( I<-I0/as.numeric( mean(diag(  as.matrix(wall) %*% K %*% t( as.matrix(wall)) )) - colMeans(wall)%*%K%*%colMeans(wall) ) )
  (I <-
      I0 / as.numeric(sum(diag(
        wall_wbar %*% K %*% t(wall_wbar)
      ))))
  
  
  return(I)
  
  
  
}




drfCI <- function(X, Y, B, sampling = "binomial",...) {
  
  n <- dim(X)[1]
  
  # compute point estimator and DRF per halfsample
  # weightsb: B times n matrix of weights
  DRFlist <- lapply(seq_len(B), function(b) {
    
    # half-sample index
    indexb <- if (sampling == "binomial") {
      seq_len(n)[as.logical(rbinom(n, size = 1, prob = 0.5))]
    } else {
      sample(seq_len(n), floor(n / 2), replace = FALSE)
    }
    
    ## Using normal Bootstrap on the data and refitting DRF
    DRFb <- 
      drf(X = X[indexb, , drop = F], Y = Y[indexb, , drop = F],
          ci.group.size = 1, ...)
    
    
    return(list(DRF = DRFb, indices = indexb))
  })
  
  return(list(DRFlist = DRFlist, X = X, Y = Y) )
}


predictdrf<- function(DRF, x, functional = NULL, ...) {
  
  ntest <- nrow(x)
  n <- nrow(DRF$Y)
  
  
  weightsb <- lapply(DRF$DRFlist, function(l) {
    
    weightsbfinal <- Matrix(0, nrow = ntest, ncol = n , sparse = TRUE)
    
    weightsbfinal[, l$indices] <- predict(l$DRF, x)$weights 
    
    return(weightsbfinal)
  })
  
  
  weightsall <- Reduce("+", weightsb) / length(weightsb)
  
  if (!is.null(functional)) {
    ## Achtung: This so far works only for one x!  
    stopifnot("Not yet implemented for several x" = ntest == 1)
    
    thetahatb <- 
      lapply(weightsb, function(w)  
        functional(weights = w , X = DRF$X, Y = DRF$Y, x = x))
    thetahatbforvar <- do.call(rbind, thetahatb)
    thetahat <- functional(weights = weightsall , X = DRF$X, Y = DRF$Y, x = x)
    thetahat <- matrix(thetahat, nrow = dim(x)[1])
    var_est <- if (dim(thetahat)[2] > 1){  
      a <- sweep(thetahatbforvar, 2, thetahat, FUN = "-")
      crossprod(a, a) / B
    } else {
      mean((c(thetahatbforvar) - c(thetahat)) ^ 2) 
    }
    
    return(list(weights = weightsall, thetahat = thetahat, weightsb = weightsb, 
                var_est = var_est ))
    
  } else {
    return(list(weights = weightsall, weightsb = weightsb ))
  }
}



condindtest <- function(X,Y,Z, x, B,...){
  ## Tests whether Y ind of Z | X=x
  
  listofargs<-list(...)
  
  if ("response.scaling" %in% names(listofargs)){
    response.scaling <- unname(unlist(listofargs["response.scaling"]))
  }else{
    # Default in DRF
    response.scaling<-T
  }
  
  
  # Fit 3 DRFs, one for (Y,Z) jointly given X, with CIs, one
  # for Y|X only and one for Z|X only 
  DRF<-drfCI(X=X,Y=cbind(Y,Z), B=B,...)
  DRF1 <- drf(X=X,Y=Y, ...)
  DRF2 <- drf(X=X,Y=Z, ...)
  
  ## predict to obtain mu(x) and mu_1(x)*mu_2(x)
  DRFpred<-predictdrf(DRF, x=x)
  w<-DRFpred$weights # mu(x)
  w1<-predict(DRF1,newdata=x)$weights #mu_1(x)
  w2<-predict(DRF2,newdata=x)$weights #mu_2(x)
  
  ## Contine here!!!
  if (response.scaling==T){
    Y.transformed <- scale(Y)
    Z.transformed <- scale(Z)
  }else{
    Y.transformed<-Y
    Z.transformed<-Z
  }
  
  bandwidth.Y <- drf:::medianHeuristic(Y.transformed)
  bandwidth.Z <- drf:::medianHeuristic(Z.transformed)
  
  # bandwidth.Y <- 1
  # bandwidth.Z <- 1
  
  k_Y<-rbfdot(sigma=bandwidth.Y)
  k_Z<-rbfdot(sigma=bandwidth.Z)
  K=kernelMatrix(k_Y, Y, y = Y)
  L=kernelMatrix(k_Z, Z, y = Z)
  
  
  ## Continue here: Idea ||m_{n,1} x \mu_{n,2} - mu_n||^2 should have the same asymptotic distribution as 
  ## ||mu_n - mu||^2 which in turn should have the same (asymptotic) distribution as ||mu_n^{S} - mu_n||^2
  
  
  # ||m_{n,1} x \mu_{n,2} - mu_n||^2 
  teststat<-diag(w%*%(K*L)%*%t(w)) + diag(w1%*%K%*%t(w1))*diag(w2%*%L%*%t(w2)) - 2* diag(w%*%((K%*%t(w1))*(L%*%t(w2))))
  # Using same DRF
  #teststat<-diag(w%*%(K*L)%*%t(w)) + diag(w%*%K%*%t(w))*diag(w%*%L%*%t(w)) - 2* diag(w%*%((K%*%t(w))*(L%*%t(w))))
  
  
  
  teststat<-matrix(teststat, ncol=nrow(x), byrow = T)
  
  nulldist<-lapply(DRFpred$weightsb, function(wb){
    
    # w: ntest x n matrix of weights
    #diag(t(w)%*%K%*%L%*%w) + diag(t(w1)%*%K%*%w1)*diag(t(w2)%*%K%*%w2) - 2* t(w)%*%((K%*%w1)*(K%*%w2))
    # ||mu_n^{S} - mu_n||^2 = || (wb-w)'%*% (k_1, ...., k_n) ||^2
    #diag(t(wb)%*%K%*%L%*%wb) + diag(t(w)%*%K%*%L%*%w) - 2* t(w)%*%((K%*%w1)*(K%*%w2))
    # ||mu_n^{S} - mu_n||^2
    diag( (wb-w)%*%(K*L)%*%t(wb-w) )
    
    
  })
  
  nulldist <- do.call(rbind,nulldist)
  
  pval <- mean(sapply(seq_len(ncol(nulldist)), function(j) teststat[,j] <= nulldist[,j] )) 
  
  #pval<-mean(teststat < nulldist)
  
  
  return(list(pval = pval, nulldist = nulldist, teststat = teststat))
  
}