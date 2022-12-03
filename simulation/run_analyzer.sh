#!/bin/bash
#SBATCH -A b1139
#SBATCH -p b1139
#SBATCH -t 04:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80G
#SBATCH --job-name="analyze_generic_PMCmode_EIR_vaccSP_IIV"
#SBATCH --error=log/slurm_%A_%a.err
#SBATCH --output=log/slurm_%A_%a.out


module purge all
module load python/anaconda3.6
source activate /home/mrm9534/environments/dtk-tools-p36
#source activate /home/mrm9534/environments/dtk-tools-p36_July

#cd /home/mrm9534/dtk-tools-p36/helper_tools/
#python local_db_fixer.py  -id 2022_07_23_11_20_33_996774 --status Succeeded
cd /home/mrm9534/gitrepos/ipti_pmc/simulations/analyzer/ 
#python analyzer_cohort.py --exp_name generic_PMC_CM_vaccSP_IIV --exp_id 2022_07_23_11_20_33_996774 --expdirsub '_pmc_rtss_generic'
python analyzer_cohort_generic2.py --exp_name mrm9534_generic_PMC3mode_shifted_vaccSP_IIV --exp_id da0caf99-8533-ed11-a9fc-b88303911bc1 --expdirsub '_pmc_rtss_generic'
#python analyzer_cohort_generic.py --exp_name mrm9534_generic_PMC3mode_shifted_vaccSP_IIV --exp_id da0caf99-8533-ed11-a9fc-b88303911bc1 --expdirsub '_pmc_rtss_generic'