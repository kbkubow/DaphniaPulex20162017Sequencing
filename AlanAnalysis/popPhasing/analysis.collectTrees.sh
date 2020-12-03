#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-02:00:00 # Running time of 4 days
#SBATCH --mem 20G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/maketree.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/maketree.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.collectTrees.sh
### sacct -u aob2x -j 19149417
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/maketree.19148067

### sbatch --array=234 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.makeTrees.sh
### sacct -u aob2x -j 19143682
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.19110025_10.err

### modules
  module load samtools parallel
  module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3

  Rscript --vanilla /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.collect.trees.R
