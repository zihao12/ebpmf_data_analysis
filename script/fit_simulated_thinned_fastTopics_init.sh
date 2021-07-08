#!/bin/bash

# The shell commands below will submit Slurm jobs to perform the
# Poisson NMF model fitting for all data sets, and for different
# choices of the model parameters and optimization settings.
SCRIPT_FIT=fit_simulated_thinned_fastTopics_init.sbatch

## on sla_simulated data
datafile=../data/sim/sla_simulated.Rds
initfile_tmp=../output/fastTopics_fit/fit_sla_fastTopics
outputfile_tmp=../output/sim_thinned/fit_sla_simualted_thinned_fastTopics_init_truth
for i in {2..20}
do
  inputfile=${initfile_tmp}_k${i}.Rds
  outputfile=${outputfile_tmp}_k${i}.Rds
  sbatch ${SCRIPT_FIT} ${datafile} ${inputfile} ${outputfile}
done

## on sla_simulated_0.20 data
datafile=../data/sim/sla_simulated_0.20.Rds
initfile_tmp=../output/fastTopics_fit/fit_sla_fastTopics
outputfile_tmp=../output/sim_thinned/fit_sla_simualted_thinned_0.20_fastTopics_init_truth
for i in {2..20}
do
  inputfile=${initfile_tmp}_k${i}.Rds
  outputfile=${outputfile_tmp}_k${i}.Rds
  sbatch ${SCRIPT_FIT} ${datafile} ${inputfile} ${outputfile}
done

## on sla_simulated_0.10 data
datafile=../data/sim/sla_simulated_0.10.Rds
initfile_tmp=../output/fastTopics_fit/fit_sla_fastTopics
outputfile_tmp=../output/sim_thinned/fit_sla_simualted_thinned_0.10_fastTopics_init_truth
for i in {2..20}
do
  inputfile=${initfile_tmp}_k${i}.Rds
  outputfile=${outputfile_tmp}_k${i}.Rds
  sbatch ${SCRIPT_FIT} ${datafile} ${inputfile} ${outputfile}
done

