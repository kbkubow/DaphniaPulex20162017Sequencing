#!/usr/bin/env bash
#
#SBATCH -J trioPhase_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 1:00:00 # Running time of 1 hours
#SBATCH --mem 12G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load intel/18.0 intelmpi/18.0 R/3.6.3
module load mathematica/11.1.1
module load parallel

#SLURM_ARRAY_TASK_ID=2

wd="/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing"
datadir="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase"
chr=$( grep "^${SLURM_ARRAY_TASK_ID}" ${datadir}/chrs.csv | cut -f2 -d',' )
echo $chr

Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbitRedo/formatData.R ${chr}
