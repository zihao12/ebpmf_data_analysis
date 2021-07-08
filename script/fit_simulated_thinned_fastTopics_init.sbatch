#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=12:00:00

# Get the command-line arguments.
datafile=${1}
initfile=${2}
outputfile=${3}
echo datafile=${datafile}
echo initfile=${initfile}
echo outputfile=${outputfile}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.6.1

# Run the R script.
echo "Running fit_simulated_fastTopics_init"
export OPENBLAS_NUM_THREADS=1
Rscript --no-save --no-restore --verbose ./fit_simulated_thinned_fastTopics_init.R ${datafile} ${initfile} ${outputfile} > ${outputfile}.out
