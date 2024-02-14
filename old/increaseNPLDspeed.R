


# Simulate 1000 observations from Y|X=x
Yx<-Y[sample(1:nrow(Y),size=1000, replace=T,prob=w),]

H=diag(apply( Yx, 2, drf:::medianHeuristic ))

#bandwidth <- drf:::medianHeuristic(Yx)
#densityvaly <- kde(Yx, eval.points = y, h=bandwidth)$estimate
densityvaly <- kde(Yx, eval.points = y, H=H)$estimate


datax<-lapply(1:nrow(weights), function(j) Y[sample(1:nrow(Y),size=1000, replace=T,prob=weights[j,]),] )
Hx<-sapply(1:nrow(weights), function(j) 1/sqrt(apply( datax[[j]], 2, drf:::medianHeuristic )) )

datax.transformed<-