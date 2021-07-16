#!/bin/bash
SCRIPT_FIT=fit_simulated_multinom.sbatch

## fit on data simulated from sla
modelname_tmp=../data/sim_multinom/sla_simulated_fastTopics

### method fastTopics
method=fastTopics
Rfile=fit_simulated_multinom_${method}.R
outname_tmp=../output/sim_multinom/fit_${method}_sla_simulated_fastTopics
for prop in 1 0.7 0.5 0.3 0.1
do
  for i in {2..20}
  do
    modelname=${modelname_tmp}_k${i}_${prop}.Rds
    outname=${outname_tmp}_k${i}_${prop}.Rds
    sbatch ${SCRIPT_FIT} ${Rfile} ${modelname} ${outname}
  done
done


### method ebpmf_wbg2
method=ebpmf_wbg2
Rfile=fit_simulated_multinom_${method}.R
outname_tmp=../output/sim_multinom/fit_${method}_sla_simulated_fastTopics
for prop in 1 0.7 0.5 0.3 0.1
do
  for i in {2..20}
  do
    modelname=${modelname_tmp}_k${i}_${prop}.Rds
    outname=${outname_tmp}_k${i}_${prop}.Rds
    sbatch ${SCRIPT_FIT} ${Rfile} ${modelname} ${outname}
  done
done

### method ebpmf_wbg2_pg
method=ebpmf_wbg2
Rfile=fit_simulated_multinom_${method}.R
outname_tmp=../output/sim_multinom/fit_${method}_sla_simulated_fastTopics
for prop in 1 0.7 0.5 0.3 0.1
do
  for i in {2..20}
  do
    modelname=${modelname_tmp}_k${i}_${prop}.Rds
    outname=${outname_tmp}_k${i}_${prop}.Rds
    sbatch ${SCRIPT_FIT} ${Rfile} ${modelname} ${outname}
  done
done


