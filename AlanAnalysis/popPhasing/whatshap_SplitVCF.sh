#!/usr/bin/env bash
#
#
#SBATCH -J popPhasing_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 2:00:00 # Running time of 2 hours
#SBATCH --mem 8G # Memory request of 12 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run with: sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_SplitVCF.sh
## sacct -u aob2x -j 9899403

#ijob -c1 -p standard -A berglandlab
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools
export PATH=$HOME/.local/bin:$PATH

### run whatshap

chr=Scaffold_1863_HRSCAF_2081
sample=April17_2018_D8_Male1

whatshap \
phase \
-o  /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.${chr}.phase.vcf \
--chromosome ${chr} \
--sample ${sample} \
/scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.${chr}.vcf.gz \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/${sample}_finalmap_mdup.bam
