#!/usr/bin/env bash
#
#
#SBATCH -J singleLocus_model # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:15 # Running time of 15 minutes
#SBATCH --mem 2G # Memory request of 2 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/singleLocus_model_slim_%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/singleLocus_model_slim_%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run command
module load gcc/7.1.0 slim

N=$( cat /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/slimuations/singleLocus_model.paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
CR=$( cat /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/slimuations/singleLocus_model.paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
AMR=$( cat /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/slimuations/singleLocus_model.paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
BMR=$( cat /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/slimuations/singleLocus_model.paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )

slim -d Ni=$N -d CRi=$CR -d AMRi=$AMR -d BMRi=$BMR /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/slimuations/singleLocus_model.slim
