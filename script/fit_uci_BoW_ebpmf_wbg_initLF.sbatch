#!/bin/bash

#SBATCH --partition=mstephens
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32G
#SBATCH --time=36:00:00

module load R/3.6.1

Rscript --no-save --no-restore --verbose fit_uci_BoW_ebpmf_wbg_initLF.R kos 20 5000 500 > ../output/uci_BoW/v0.4.2/fit_kos_ebpmf_wbg_initLF_K20_maxiter_5000.Rout 2>&1
