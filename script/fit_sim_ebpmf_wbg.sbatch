#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=36:00:00

module load R/3.6.1

Rscript --no-save --no-restore --verbose fit_sim_ebpmf_wbg.R sim_bg_block_n1100_p2100_K50 init.sim_bg_block_n1100_p2100_K50.Rds 5000 500 > ../output/sim/v0.4.5/fit_sim_bg_block_n1100_p2100_K50_ebpmf_wbg_maxiter_5000.Rout 2>&1
