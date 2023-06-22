
evalsynthetic <- function(dataset, L=10, n, p, num.trees, ... ){
  
  resmat <- matrix(NaN, nrow=p, ncol=L)
  rownames(resmat) <- paste0("X",1:p)
  
  for (l in 1:L){
    
    cat(l)
    
    tmp<-genData(dataset = dataset, n = n, p = p)
    
    
    X<-tmp$X
    Y<-as.matrix(tmp$y)
    colnames(X) <- paste0("X",1:ncol(X))
    
    
    
    ressynth<-drfwithVI(X, Y, B=1, num.trees=num.trees)
    
    resmat[,l] <- ressynth$VI
  }
  
  res <- rowMeans(resmat)
  resordered <- round(sort(res, decreasing = T),3)
  nams<-names(resordered)
  resordered <- matrix(resordered,nrow=1)
  colnames(resordered) <- nams
  
  
  ## Create latex table
  print(resordered %>%
    kbl(caption=paste0("Ordering for the ", dataset, " dataset") ,
        format="latex",
        #col.names = c("$I_n^{(-j)}$"),
        col.names = nams,
        align="r") %>%
    kable_minimal(full_width = F,  html_font = "Source Sans Pro"))
  
}



evalrealdata <- function(dataset, n, p = 10, ntest=round(1/3*n),... ){
  
  
  
  tmp<-genData(dataset = dataset, n = n, p = 10, ...)
  
  
  X<-tmp$X
  Y<-as.matrix(tmp$y)
  colnames(X) <- paste0("X",1:ncol(X))
  
  Xtest <- X[(round(n - ntest) + 1):n, , drop = F]
  Ytest <- Y[(round(n - ntest) + 1):n, , drop = F]
  #
  X <- X[1:round(n - ntest), , drop = F]
  Y <- Y[1:round(n - ntest), , drop = F]
  
  
  
  whichvar<-rep(NA,ncol(X)-1)
  eval<-rep(NA,ncol(X)-1)
  p<-ncol(X)
  
  for (j in 1:(p - 1)  ){
    
    
    # remove variable with smallest VI
    resreal<-featureeliminnation(X,Y)
    X<-resreal$Xnew
    Xtest<-Xtest[,!( colnames(Xtest) %in% resreal$which), drop=F]
    
    
    whichvar[j] <-resreal$which
    eval[j]<-distpredicteval(X,Y,Xtest, Ytest,dNPLD, num.trees=num.trees)
  }
  
  
  return(list(whichvar=whichvar, eval=eval))
  
}