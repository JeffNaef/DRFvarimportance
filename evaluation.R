
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