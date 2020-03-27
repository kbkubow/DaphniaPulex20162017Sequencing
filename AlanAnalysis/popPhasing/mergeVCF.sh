#!/usr/bin/env bash
#
#
#SBATCH -J popPhasing_mergeVCF # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0:30:00 # Running time of 15 minutes
#SBATCH --mem 1G # Memory request of 4GGB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# ijob -c1 -p standard -A berglandlab
### run with: sbatch --array=1 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/mergeVCF.sh
### sacct -u aob2x -j 10016822
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_whatshapp.10007470_20.out

### load modules
  module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822

## get job
  # SLURM_ARRAY_TASK_ID=1
  chr=$( cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' | sort | uniq | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )

## make sure that all vcf files are there, should be 511 pulex. 509 worked because of some diferences in naming conventions

## bgzip vcf files
#for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/*.${chr}.phase.vcf; do
#  echo "File -> $f"
#  bgzip \
#  -c \
#  -@ 1 \
#  ${f} > ${f}.gz
#done

#for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/*.${chr}.phase.vcf.gz; do
#  echo "File -> $f"
#  tabix \
#  -p vcf \
#  ${f}
#done

  ls /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/*.${chr}.phase.vcf.gz > /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${chr}.list

  bcftools \
  merge \
  -l /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${chr}.list \
  -o  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/${chr}.whatshapp.vcf \
  -O v \
  --threads 10
