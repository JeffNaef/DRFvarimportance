
evalsynthetic <- function(dataset, L=10, B=1, n, p, num.trees, MRF=F, ... ){
  
  
  if (dataset=="motivatingexample"){
    p<-2
  }
  
  if (dataset=="GP"){
    num.features<-50
  }else{
    num.features <- 10
  }
  
  resmat <- matrix(NaN, nrow=p, ncol=L)
  rownames(resmat) <- paste0("X",1:p)
  
  if (MRF==T){
    
    resmatMRF1 <- matrix(NaN, nrow=p, ncol=L)
    rownames(resmatMRF1) <- paste0("X",1:p)
    
    
    resmatMRF2 <- matrix(NaN, nrow=p, ncol=L)
    rownames(resmatMRF2) <- paste0("X",1:p)
  }
  
  
  for (l in 1:L){
    
    cat(l)
    
    tmp<-genData(dataset = dataset, n = n, p = p)
    
    
    X<-tmp$X
    Y<-as.matrix(tmp$y)
    colnames(X) <- paste0("X",1:ncol(X))
    
    if (dataset=="GP" & l==1){
      
      plot(seq(-5, 5, length.out = length(Y[1,])),Y[1,], type="l", cex=0.8, col="darkblue", ylim=c(-3,3), xlab="t", ylab="Y")
      
      for (i in 2:10){
        
        lines(seq(-5, 5, length.out = length(Y[1,])),Y[i,], cex=0.8, col="darkblue")
      }
      
      
    }  
    
    
    
    ressynth<-drfwithVI(X, Y, B=B, num.trees=num.trees, num.features=num.features, ...)
    
    resmat[,l] <- ressynth$VI
    
    if (MRF==T){
      
      MRFVI_1<-MeanOutcomeDifference(X, Y, num_trees=num.trees, ...)
      MRFVI_2<-MeanSplitImprovement(X, Y, num_trees=num.trees, ...)
      
      
      resmatMRF1[,l] <- MRFVI_1$VI
      resmatMRF2[,l] <- MRFVI_2$VI
    }
    
    
    
  }
  
  res <- rowMeans(resmat)
  resordered <- round(sort(res, decreasing = T),3)
  nams<-names(resordered)
  resordered <- matrix(resordered,nrow=1)
  colnames(resordered) <- nams
  
  
  if (MRF==T){
    resMRF1 <- rowMeans(resmatMRF1)
    resorderedMRF1 <- round(sort(resMRF1, decreasing = T),3)
    nams<-names(resorderedMRF1)
    resorderedMRF1 <- matrix(resorderedMRF1,nrow=1)
    colnames(resorderedMRF1) <- nams
    
    resMRF2 <- rowMeans(resmatMRF2)
    resorderedMRF2 <- round(sort(resMRF2, decreasing = T),3)
    nams<-names(resorderedMRF2)
    resorderedMRF2 <- matrix(resorderedMRF2,nrow=1)
    colnames(resorderedMRF2) <- nams
    
  }
  
  
  ## Create latex table
  print(resordered %>%
          kbl(caption=paste0("Ordering for the ", dataset, " dataset") ,
              format="latex",
              #col.names = c("$I_n^{(-j)}$"),
              col.names = nams,
              align="r") %>%
          kable_minimal(full_width = F,  html_font = "Source Sans Pro"))
  
}



evalrealdata <- function(dataset, n, B=1, ntest=round(1/2*n), num.trees=500, successivVI=F , ... ){
  
  
  
  tmp<-genData(dataset = dataset, n = n, ...)
  
  
  X<-tmp$X
  Y<-as.matrix(tmp$y)
  
  if (is.null(colnames(X))){
    colnames(X) <- paste0("X",1:ncol(X))
  }
  
  
  Xtest <- X[(round(n - ntest) + 1):n, , drop = F]
  Ytest <- Y[(round(n - ntest) + 1):n, , drop = F]
  #
  X <- X[1:round(n - ntest), , drop = F]
  Y <- Y[1:round(n - ntest), , drop = F]
  
  
  
  
  
  if (successivVI==T){
    
    whichvar<-rep(NA,ncol(X)-1)
    eval<-rep(NA,ncol(X)-1)
    p<-ncol(X)
    
    
    ##Calculate variable importance several times
    for (j in 1:(p - 1)  ){
      
      
      # remove variable with smallest VI
      resreal<-featureeliminnation(X,Y)
      X<-resreal$Xnew
      Xtest<-Xtest[,!( colnames(Xtest) %in% resreal$which), drop=F]
      
      
      whichvar[j] <-resreal$which
      eval[j]<-distpredicteval(X,Y,Xtest, Ytest,dNPLD, num.trees=num.trees)
      #eval[j]<-distpredicteval(X,Y,Xtest, Ytest,dMMD, num.trees=num.trees)
    }
    
    
    return(list(whichvar=whichvar, eval=eval))
    
  }else{
    
    
    
    ##Calculate variable importance only once
    ressynth<-drfwithVI(X, Y, B=B, num.trees=num.trees, ...)
    
    
    return( list(ressynth=ressynth, X=X, Y=Y, Xtest=Xtest, Ytest=Ytest)   )
  }
  
}




evalfinal <- function(Vimp, X ,Y ,Xtest, Ytest, metrics=c("MMD","NPLD","MAD"), num.trees ){
  
  if (ncol(Ytest) > 1 & "MAD" %in% metrics){
    metrics <- metrics[!( metrics %in% "MAD") ]
  }
  
  Vimp<-sort(Vimp)
  
  if ( is.null(names(Vimp)) ){
    stop("Need names for later")  
  }
  
  
  evalMMD<-matrix(0, nrow=ncol(X))
  evalNPLD<-matrix(0, nrow=ncol(X))
  evalMAD<-matrix(0, nrow=ncol(X))
  
  ###Idea: Create a function that takes a variable importance measure and does this loop!!
  
  for (j in 1:ncol(X)){
    
    print(paste0('Running importance for variable X', j, '...'))
    
    
    if (j==1){
      
      if ("MMD" %in% metrics){
        
        evalMMD[j]<-distpredicteval(  X,Y,
                                      Xtest, 
                                      Ytest,d="MMD", num.trees=num.trees)
      }
      
      if ("NPLD" %in% metrics){
        evalNPLD[j]<-distpredicteval(  X,Y,
                                       Xtest, 
                                       Ytest,d="NPLD", num.trees=num.trees) 
      }
      
      
      DRFall <- drf(X=X, Y=Y, num.trees=num.trees)
      quantpredictall<-predict(DRFall, newdata=Xtest, functional="quantile",quantiles=c(0.5))
      evalMAD[j] <- mean(sapply(1:nrow(Xtest), function(j)  abs(Ytest[j] - quantpredictall$quantile[,,"q=0.5"][j]) ))
      
      
    }else{
      
      
      if ("MMD" %in% metrics){
        evalMMD[j]<-distpredicteval(  X[,!(colnames(X) %in% names(Vimp[1:(j-1)])), drop=F],Y,
                                      Xtest[,!(colnames(Xtest) %in% names(Vimp[1:(j-1)])), drop=F], 
                                      Ytest,d="MMD", num.trees=num.trees)
      }
      
      if ("NPLD" %in% metrics){
        evalNPLD[j]<-distpredicteval(  X[,!(colnames(X) %in% names(Vimp[1:(j-1)])), drop=F],Y,
                                       Xtest[,!(colnames(Xtest) %in% names(Vimp[1:(j-1)])), drop=F], 
                                       Ytest,d="NPLD", num.trees=num.trees)
      }
      
      
      if ("MAD" %in% metrics){
        DRFall <- drf(X=X[,!(colnames(X) %in% names(Vimp[1:(j-1)])), drop=F], Y=Y, num.trees=num.trees)
        quantpredictall<-predict(DRFall, newdata=Xtest[,!(colnames(Xtest) %in% names(Vimp[1:(j-1)])), drop=F], functional="quantile",quantiles=c(0.5))
        evalMAD[j] <- mean(sapply(1:nrow(Xtest), function(j)  abs(Ytest[j] - quantpredictall$quantile[,,"q=0.5"][j]) ))
      }
      
    }
    
  }
  
  return(list(Vimp=Vimp, evalMMD=evalMMD, evalNPLD=evalNPLD, evalMAD=evalMAD ))
  
}




# Compute the MMD
mmd <- function(X, Y, sigma) {
  # Define the Gaussian kernel
  k_X <- rbfdot(sigma = sigma)
  k_Y <- rbfdot(sigma = sigma)
  
  # Compute the kernel matrices
  K_X <- kernelMatrix(k_X, X, y = X)
  K_Y <- kernelMatrix(k_Y, Y, y = Y)
  K_XY <- kernelMatrix(k_X, t(X), y = Y)
  
  # Compute the average kernel values
  XX <- mean(K_X)
  YY <- mean(K_Y)
  XY <- mean(K_XY)
  
  # Compute the MMD
  sqrt(XX + YY - 2 * XY)
}



## Define a distance function D
dMMD<-function(w, Y,y, K_Y){
  
  # Simulate 1000 observations from Y|X=x
  #Yx<-Y[sample(1:nrow(Y),size=1000, replace=T,prob=w),]
  
  
  sigma <- drf:::medianHeuristic(Y)
  ## Can do this better by direct calculation!
  
  
  # Define the Gaussian kernel
  k_y <- rbfdot(sigma = sigma)
  k_Y <- rbfdot(sigma = sigma)
  
  # Compute the kernel matrices
  K_y <- kernelMatrix(k_y, y, y = y)
  #K_Y <- kernelMatrix(k_Y, Y, y = Y)
  K_yY <- kernelMatrix(k_y, t(y), y = Y)
  
  # Compute the average kernel values
  XX <- K_y
  YY <- w%*%K_Y%*%w
  XY <- sum(w*K_yY)
  
  # Compute the MMD
  
  
  return(sqrt(XX + YY - 2 * XY))
  
}

dNPLD <- function(w,Y,y){
  
  
  # Simulate 1000 observations from Y|X=x
  Yx<-Y[sample(1:nrow(Y),size=1000, replace=T,prob=w),]
  
  #bandwidth <- drf:::medianHeuristic(Yx)
  #densityvaly <- kde(Yx, eval.points = y, h=bandwidth)$estimate
  densityvaly <- kde(Yx, eval.points = y)$estimate
  
  return( - log(densityvaly))
  
}





distpredicteval <- function(X,Y,Xtest, Ytest,d="MMD", parallel=F, ...){
  
  
  ## d: a function that takes in d, Y and a test point y and evaluates a metric
  ## of accuracy
  
  if (ncol(Xtest)!=ncol(X)){
    stop("Xtest must have same number of columns as X")
  }
  
  
  # Step 1: Fit and Predict
  DRF<-drf(X,Y, ...)
  weights<-predict(DRF,newdata=Xtest)$weights
  
  
  sigmatest <- drf:::medianHeuristic(Ytest)
  sigmatrain <- drf:::medianHeuristic(Y)
  
  if (d=="MMD"){
    
    Y.transformed <- scale(Y)
    Ytest.transformed <- scale(Ytest)
    
    # Define the Gaussian kernel
    k_y <- rbfdot(sigma = sigmatest)
    k_Y <- rbfdot(sigma = sigmatrain)
    
    # Compute the kernel matrices
    K_y <- kernelMatrix(k_y, Ytest.transformed, y = Ytest.transformed)
    K_Y <- kernelMatrix(k_Y, Y.transformed, y = Y.transformed)
    K_yY <- kernelMatrix(k_y, Ytest.transformed, y = Y.transformed)
    
    
    
    yy<-diag(K_y)
    YY<-diag(weights%*%K_Y%*%t(weights))
    Yy<-diag(K_yY%*%t(weights))
    
    res <- mean(sqrt(yy+YY-2*Yy))
    
  }else if (d=="NPLD"){
    
    if (parallel==T){
      
      D <- foreach(i = 1:nrow(weights), .export=c("dNPLD","mmd"), .packages = c("ks", "kernlab", "drf", "Matrix"), .combine=rbind) %dopar% {
        
        result <- - dNPLD(weights[i,], Y, Ytest[i,])  # Make sure you replace 'your_package_containing_d_function' with the actual package name
        return(result)
      }
    }else{
      
      D<-sapply(1:nrow(weights),function(i){
        
        if(i%%10==0){cat(i)}
        
        dNPLD(weights[i,], Y, Ytest[i,])
        
      } )
    }
    
    res<-mean(D, na.rm=T, trim=0.05)
    
  }
  
  
  
  return(res)
  
}

