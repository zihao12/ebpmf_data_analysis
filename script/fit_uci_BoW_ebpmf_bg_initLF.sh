#!/bin/bash

Rscript --no-save --no-restore --verbose fit_uci_BoW_ebpmf_bg_initLF.R kos 20 5000 500 > ../output/uci_BoW/v0.3.9/fit_kos_ebpmf_bg_initLF_K20_maxiter_5000.Rout 2>&1
