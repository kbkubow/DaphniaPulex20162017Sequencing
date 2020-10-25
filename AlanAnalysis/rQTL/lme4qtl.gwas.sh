#!/usr/bin/env bash
#
#SBATCH -J lme4qtl # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=10 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 7:00:00 # Running time of 1 hours
#SBATCH --mem 4G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run as
# sbatch --array=1 ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/lme4qtl.gwas.sh
# sacct -j 17584183
# cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.17584183_9.out

module load intel/18.0 intelmpi/18.0 R/3.6.3

#SLURM_ARRAY_TASK_ID=1
wd="/scratch/aob2x/daphnia_hwe_sims"

Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/lme4qtl.gwas.R ${SLURM_ARRAY_TASK_ID}
