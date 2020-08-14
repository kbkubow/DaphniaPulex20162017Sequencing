#!/usr/bin/env bash
#
#SBATCH -J split_and_run # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2:00:00 ### most jobs should run in 60 minutes or less; the mitochondria takes a lot longer to run through pool-snp
#SBATCH --mem 5G
##SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
##SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/treeMix/runTreeMix.sh
# sacct -j 14267973

module load intel/18.0  intelmpi/18.0 gsl/2.4 R/3.6.3 boost/1.60.0; R


~/treemix/src/treemix \
-i /scratch/aob2x/treemixIn.pond_species.gz \
-k 500 \
-m 2 \
-o /scratch/aob2x/treemixIn.pond_species \
-root pulicaria_Pond21
