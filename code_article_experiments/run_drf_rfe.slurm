#!/bin/bash

set -e
DATASETS=("enb" "jura" "wq" "scm20d" "rf1" 
          "oes97" "air" "births2" "wage")
for DATASET in ${DATASETS[@]}; do
	export DATASET=$DATASET
	sbatch submit_rfe.slurm
done