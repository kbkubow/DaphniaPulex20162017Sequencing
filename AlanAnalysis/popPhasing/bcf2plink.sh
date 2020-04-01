#!/usr/bin/env bash

#SBATCH -J bcf2plink # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 00:10:00 # Running time of 4 days
#SBATCH --mem 20G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/bcf2plink.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/bcf2plink.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run with: sbatch --array=1-12 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/bcf2plink.sh
### sacct -u aob2x -j 10127549
#
## ijob -c1 -p standard -A berglandlab
#SLURM_ARRAY_TASK_ID=1

chr=$( cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' | sort | uniq | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )

module load plink/1.90b6.16

plink \
--bcf /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/${chr}.whatshapp.onePerSC.bcf \
--make-bed \
--double-id \
--allow-extra-chr \
--out /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/${chr}.whatshapp.onePerSC
