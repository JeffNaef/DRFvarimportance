Instructions to run experiments

1. Simulated data
Experiments with simulated data can be reproduced using the R scripts 'drf_vimp_simulated_xp_univariate.R',
and 'drf_vimp_simulated_xp_bivariate_functional.R' (at the beginning of the script, set 'xp_num <- 2' for the bivariate case, 
and 'xp_num <- 3' for the functional case). See required R packages at the beginning of each script.

2. RFE for real datasets
Experiments can be run on a cluster with a slurm system by running the batch file 'run_drf_rfe.slurm', 
which launches a separate slurm job for each dataset using 'submit_rfe.slurm'.
Then, the R script 'drf_rfe.R' is executed, and results are stored in a Rdata file for each dataset in the folder 'mtr_results'.
Next, the script 'plot_mtr_results.R' can be used to generate plots, stored in the folder 'mtr_plots'.
See required R packages at the beginning of each R script.