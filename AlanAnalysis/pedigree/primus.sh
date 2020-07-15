#!/usr/bin/env bash
#
#
#SBATCH -J primus # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 6:00:00 # Running time of 15 minutes
#SBATCH --mem 9G # Memory request of 4GGB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/primus.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/primus.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run primus: sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/pedigree/primus.sh
### sacct -j 10077999

### convert vcf file
module load plink/1.90b6.16 bcftools/1.9


### rename chrs
bcftools \
annotate \
--rename-chrs /scratch/aob2x/daphnia_hwe_sims/pedigree/chr_rename.delim \
-O b \
-o /scratch/aob2x/daphnia_hwe_sims/pedigree/clean/onePerSC.LDprune.bcf \
/scratch/aob2x/daphnia_hwe_sims/pedigree/clean/onePerSC.LDprune.vcf



plink \
--bcf /scratch/aob2x/daphnia_hwe_sims/pedigree/clean/onePerSC.LDprune.bcf \
--make-bed \
--maf 0.001 \
--make-founders \
--double-id \
--allow-extra-chr \
--out /scratch/aob2x/daphnia_hwe_sims/pedigree/clean/onePerSC.LDprune


plink \
--bfile /scratch/aob2x/daphnia_hwe_sims/pedigree/clean/onePerSC.LDprune \
--genome \
--out /scratch/aob2x/daphnia_hwe_sims/pedigree/clean/onePerSC.LDprune

less -S /scratch/aob2x/daphnia_hwe_sims/pedigree/clean/onePerSC.LDprune.genome

wd=/scratch/aob2x/daphnia_hwe_sims/pedigree/PRIMUS_v1.9.0/bin


# rm -fr /scratch/aob2x/daphnia_hwe_sims/pedigree/primusOut/
mkdir /scratch/aob2x/daphnia_hwe_sims/pedigree/primusOut/clean/out

${wd}/run_PRIMUS.pl \
--file /scratch/aob2x/daphnia_hwe_sims/pedigree/clean/onePerSC.LDprune \
--genome \
--age_file /scratch/aob2x/daphnia_hwe_sims/pedigree/snprelatre_output.ages.delim \
-o /scratch/aob2x/daphnia_hwe_sims/pedigree/primusOut/clean/out \
-v 3








${wd}/run_PRIMUS.pl \
-i FILE=/scratch/aob2x/daphnia_hwe_sims/pedigree/clean/onePerSC.LDprune.genome \
IBD0=7 IBD1=8 IBD2=9 PI_HAT=10 \
-o /scratch/aob2x/daphnia_hwe_sims/pedigree/primusOut/clean/out \
--degree_rel_cutoff 1 \
--max_gen_gap 2 \
-v 2
