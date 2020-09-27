#!/bin/bash

module load R/3.6.1

Rscript --no-save --no-restore --verbose fit_uci_BoW_ebpmf_bg.R kos 20 5000 500 > ../output/uci_BoW/v0.3.9/fit_kos_ebpmf_bg_K20_maxiter_5000.Rout 2>&1
