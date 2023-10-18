
# load packages
library(ggplot2)

datasets_list <- c("enb", "jura", "wq", "air", "births2", "rf1", "scm20d",
                   "wage", "oes97")

for (dataset in datasets_list){
  
  # load results
  load(paste0('mtr_results/rfe_results_', dataset, '.Rdata'))
  p <- nrow(mmd_rfe_drf)
  niter <- ncol(mmd_rfe_drf)
  
  # sd plots
  mmd_rfe_drf_mean <- round(apply(mmd_rfe_drf, 1, mean), 4)
  mmd_rfe_drf_sd <- round(apply(mmd_rfe_drf, 1, sd)/sqrt(niter), 4)
  mmd_rfe_vimp_mean <- round(apply(mmd_rfe_vimp, 1, mean), 4)
  mmd_rfe_vimp_sd <- round(apply(mmd_rfe_vimp, 1, sd)/sqrt(niter), 4)
  
  df1 <- data.frame(cbind(rep("drf", p), rep(1:p), mmd_rfe_drf_mean, mmd_rfe_drf_sd))
  df2 <- data.frame(cbind(rep("vimp", p), rep(1:p), mmd_rfe_vimp_mean, mmd_rfe_vimp_sd))
  colnames(df1) <- colnames(df2) <- c('algo', 'dim', 'mmd', 'sd')
  df_all <- rbind(df1, df2)
  df_all[,2:4] <- apply(df_all[,2:4], 2, as.numeric)
  
  pdodge <- position_dodge(0.2)
  rfe_plot <- ggplot(df_all[df_all$dim >= 3 & df_all$dim <= 100,], aes(x=dim, y=mmd, colour=algo)) +
    geom_errorbar(aes(ymin=mmd-sd, ymax=mmd+sd), width=0.8, size = 0.5, position=pdodge) +
    geom_point(position=pdodge, size = 1) + 
    xlab('Number of Variables') + ylab('MMD loss') +
    theme_classic() +
    theme(legend.position = "none") +
    ggtitle(paste0('RFE for \"', dataset, '\" dataset')) +
    theme(axis.line.x = element_line(colour = 'black'),
          axis.line.y = element_line(colour = 'black'),
          text = element_text(size=14),
          plot.title = element_text(hjust = 0.5, size=18, face="italic"))
  
  rfe_plot
  ggsave(paste0('mtr_plots/RFE_mmd_sd_', dataset, '.png'), rfe_plot, width = 9, height = 6)
  
}
