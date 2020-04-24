#!/usr/bin/env bash
#
##SBATCH -J truffle_pulexOnly # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-03:00:00 # Running time of 5 days
#SBATCH --mem 100G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/truffle_pulexOnly.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/truffle_pulexOnly.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# ijob -c1 -p standard -A berglandlab
### run with: sbatch --array=1-12 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/mergeVCF.sh
### sacct -u aob2x -j 10259348
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.10096702_1.err

module load samtools htslib bcftools/1.9

### pulex only samples
  cat /project/berglandlab/Karen/MappingDec2019/CloneInfoFilePulexandObtusa_withmedrd_update20200324 | grep "pulex" | cut -f1 > \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/pulexSamples.list

### make pulex only file
  bcftools view \
  -S /scratch/aob2x/daphnia_hwe_sims/popPhase/pulexSamples.list \
  -Ob \
  --threads 20 \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.vcf.gz > \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.pulexOnly.vcf.gz

### truffle
  twd="/scratch/aob2x/daphnia_hwe_sims/popPhase/truffle"

  ${twd}/truffle \
  --vcf /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.vcf.gz \
  --maf 0.05 \
  --cpu 20  \
  --segments
