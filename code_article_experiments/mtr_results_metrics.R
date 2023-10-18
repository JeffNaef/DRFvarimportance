
library(xtable)

datasets_list <- c("enb", "jura", "wq", "air", "births2", "rf1", "scm20d",
                   "wage", "oes97")

df_metrics <- sapply(datasets_list, function(dataset){
  
  # load results
  load(paste0('mtr_resultsrfe_results_', dataset, '.Rdata'))
  p <- nrow(mmd_rfe_drf)
  niter <- ncol(mmd_rfe_drf)
  d <- ncol(data$Y)
  
  # aggregated metrics
  mmd_rfe_drf_sum <- apply(mmd_rfe_drf[1:min(100, p),], 2, sum)
  mmd_rfe_vimp_sum <- apply(mmd_rfe_vimp[1:min(100, p),], 2, sum)
  
  c(min(2000, nsample), p, d, mean(mmd_rfe_vimp_sum), sd(mmd_rfe_vimp_sum)/sqrt(niter), mean(mmd_rfe_drf_sum), sd(mmd_rfe_drf_sum)/sqrt(niter))
  
})

df_metrics <- as.data.frame(t(df_metrics))
colnames(df_metrics) <- c('n', 'p', 'd', 'vimp', 'vimp_std', 'drf', 'drf_std')
df_metrics$delta <- df_metrics$drf - df_metrics$vimp
df_metrics$std_sum <- df_metrics$drf_std + df_metrics$vimp_std
df_metrics[3:8] <- round(df_metrics[3:8], 3)
df_metrics

xtable(df_metrics[1:7], digits = c(0, 0, 0, 0, 3, 3, 3, 3))
