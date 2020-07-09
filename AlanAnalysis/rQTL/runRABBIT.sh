#!/usr/bin/env bash
#
#
#SBATCH -J trioPhase_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 12:00:00 # Running time of 2 hours
#SBATCH --mem 12G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load mathematica

math -noprompt -script /scratch/aob2x/daphnia_hwe_sims/trioPhase/AxC_Rabbit.m
