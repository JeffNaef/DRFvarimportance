#!/bin/bash

#SBATCH -p skl.q,ccl.q,icl.q
#SBATCH -c 21
#SBATCH -t 24:00:00
#SBATCH -J R-DRF
#SBATCH -o logs/log_reg_%j.o
#SBATCH -e logs/log_reg_%j.e


Rscript --vanilla drf_rfe.R $DATASET 40 FALSE