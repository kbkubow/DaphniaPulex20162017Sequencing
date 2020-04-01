#!/usr/bin/env bash
#SBATCH -J copy_bams
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-4:00:00
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/popPhase/slurmOut/copy_bams.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/popPhase/slurmOut/copy_bams.%A_%a.err # Standard error

# ijob -c1 -p standard -A berglandlab
# submit as: sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/copyBams.sh

cp -r /project/berglandlab/Karen/MappingDec2019/bams/* \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/.

cp /project/berglandlab/Karen/MappingDec2019/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf \
/scratch/aob2x/daphnia_hwe_sims/popPhase/.
