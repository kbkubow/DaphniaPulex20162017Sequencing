#!/usr/bin/env bash
#
#SBATCH -J multipool # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/multipool.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/multipool.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc/7.1.0  openmpi/3.1.4 python/2.7.16

# SLURM_ARRAY_TASK_ID=1
# run as: nJobs=$( ls /scratch/aob2x/daphnia_hwe_sims/multipool/inputData/* | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | sort | uniq | wc -l )
# run as: sbatch --array=1-${nJobs} 

chr=$( ls /scratch/aob2x/daphnia_hwe_sims/multipool/inputData/* | rev | cut -f1 -d'/' | rev | cut -f1 -d'.' | sort | uniq | \
awk -v currJob=${SLURM_ARRAY_TASK_ID} '{if(currJob==NR) print $0}' )

/scratch/aob2x/daphnia_hwe_sims/multipool/multipool/mp_inference.py -n 100 -m contrast \
-o /scratch/aob2x/daphnia_hwe_sims/multipool/output/${chr}.out \
/scratch/aob2x/daphnia_hwe_sims/multipool/inputData/${chr}.male.count.delim \
/scratch/aob2x/daphnia_hwe_sims/multipool/inputData/${chr}.pe.count.delim
