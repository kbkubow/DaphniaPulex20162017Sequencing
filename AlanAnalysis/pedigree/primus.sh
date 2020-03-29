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


module load perl

${wd}/run_PRIMUS.pl \
-i FILE=/scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim \
IBD0=7 IBD1=8 IBD2=9 PI_HAT=10 \
-o /scratch/aob2x/daphnia_hwe_sims/pedigree/primusOut/
-v 3 &
