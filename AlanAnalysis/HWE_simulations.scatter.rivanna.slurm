#!/usr/bin/env bash
# slurm_stouffers.slurm
#
#SBATCH -J hwe_simulations # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-00:90 # Running time of 30 minutes
#SBATCH --mem 8 # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/diapause_gwas/slurmOut/HWE_simulations.scatter.rivanna_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/diapause_gwas/slurmOut/HWE_simulations.scatter.rivanna%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run command
module load gcc R/3.5.1
echo Rscript /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/HWE_simulations.scatter.rivanna.R "$SLURM_ARRAY_TASK_ID"
Rscript /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/HWE_simulations.scatter.rivanna.R "$SLURM_ARRAY_TASK_ID"
