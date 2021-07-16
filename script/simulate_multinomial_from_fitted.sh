#!/bin/bash

SCRIPT_FIT=simulate_multinomial_from_fitted.R

## on sla_simulated data
modelname_tmp=../output/fastTopics_fit/fit_sla_fastTopics
outname_tmp=../data/sim_multinom/sla_simulated_fastTopics

for prop in 1 0.7 0.5 0.3 0.1
do
  for i in {2..20}
  do
    modelname=${modelname_tmp}_k${i}.Rds
    outname=${outname_tmp}_k${i}_${prop}.Rds
    Rscript ${SCRIPT_FIT} ${modelname} ${outname} ${prop}
  done
done


