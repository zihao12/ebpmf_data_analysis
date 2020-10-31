#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=36:00:00

module load R/3.6.1

Rscript --no-save --no-restore --verbose fit_sim_pmf_from_truth.R sim_bg_block_n1100_p2100_K50 5000 1000 1 > ../output/sim/v0.4.5/exper1/fit_sim_bg_block_n1100_p2100_K50_pmf_maxiter_5000_from_truth.Rout 2>&1
