#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:10:00 # Running time of 4 days
#SBATCH --mem 2G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=8,9 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/getBamSlices_Daphnia00787/RNA_bamslices_covergeDepth.sh
### sacct -u aob2x -j 21042932
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.21042932_9.out

module load samtools gcc/9.2.0  openmpi/3.1.6 python/3.7.7 bedtools/2.29.2

#SLURM_ARRAY_TASK_ID=9
wd=/scratch/aob2x/daphnia_hwe_sims/

echo ${SLURM_ARRAY_TASK_ID}
echo $samp
