#!/bin/bash
#SBATCH -A b1139
#SBATCH -p b1139
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --job-name="analyze"
#SBATCH --error=log/errors/analyze.%j.err
#SBATCH --output=log/outputs/analyze.%j.out


module purge all
ml R/4.1.1
cd /home/mrm9534/gitrepos/ipti_pmc/
Rscript pmc_rtss_generic/run_aggregate_generic.R 'generic_PMCsingle_decayshape'

