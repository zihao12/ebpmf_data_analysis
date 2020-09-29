#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=36:00:00

module load R/3.6.1

Rscript --no-save --no-restore --verbose fit_sla_ebpmf_wbg_initLF.R sla 10 5000 500 > ../output/SLA/v0.4.5/fit_sla_ebpmf_wbg_K10_maxiter_5000.Rout 2>&1
