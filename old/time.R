


library(MulvariateRandomForestVarImp)
## basic example code
set.seed(49)

n<-200

X <- matrix(runif(n*5), n, 5)
Y <- matrix(runif(n*2), n, 2)


## DRF brute force
start.time <- Sys.time()
ressynth<-drfwithVI(X, Y, B=1, num.trees=2000, num.features=10)
ressynth$VI

end.time <- Sys.time()

end.time-start.time





## MRF

start.time <- Sys.time()

split_improvement_importance <- MeanSplitImprovement(X, Y)
split_improvement_importance
#> [1] 0.8066173 2.8909635 3.4591123 0.6227943 0.5138745

end.time <- Sys.time()

end.time-start.time


start.time <- Sys.time()

mean_outccome_diff_importance <- MeanOutcomeDifference(X, Y)
mean_outccome_diff_importance
#>           [,1]      [,2]
#> [1,] 0.2458139 0.3182474
#> [2,] 0.2712269 0.2915053
#> [3,] 0.2125802 0.2023291
#> [4,] 0.2819759 0.2519035
#> [5,] 0.1238451 0.1958629

end.time <- Sys.time()

end.time-start.time