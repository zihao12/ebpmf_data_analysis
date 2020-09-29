#!/bin/bash

#SBATCH --partition=mstephens
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=36:00:00

module load R/3.6.1

Rscript --no-save --no-restore --verbose fit_uci_BoW_ebpmf_bg_initLF.R kos 20 5000 500 > ../output/uci_BoW/v0.3.9/fit_kos_ebpmf_bg_initLF_K20_maxiter_5000.Rout 2>&1
