#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=36:00:00

# Get the command-line arguments.
method=${1}
datafile=${2}
initfile=${3}
outputfile=${4}
echo method=${method}
echo datafile=${datafile}
echo initfile=${initfile}
echo outputfile=${outputfile}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.6.1

# Run the R script.
echo "Running fit_simulated_${method}_init"
export OPENBLAS_NUM_THREADS=1
Rscript --no-save --no-restore --verbose ./fit_simulated_thinned_${method}_init.R ${datafile} ${initfile} ${outputfile} > ${outputfile}.out
