#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=36:00:00

# Get the command-line arguments.
Rfile=${1}
datafile=${2}
outputfile=${3}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.6.1

# Run the R script.
export OPENBLAS_NUM_THREADS=1
echo Rscript --no-save --no-restore --verbose ${Rfile} ${datafile} ${outputfile} > ${outputfile}.out
Rscript --no-save --no-restore --verbose ${Rfile} ${datafile} ${outputfile} > ${outputfile}.out
