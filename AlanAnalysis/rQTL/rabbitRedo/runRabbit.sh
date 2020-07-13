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


### run as
# sbatch --array=1-12%1 ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbitRedo/runRabbit.sh
# sacct -j 13247208
# cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.13247110_2.out


module load intel/18.0 intelmpi/18.0 R/3.6.3
module load mathematica
module load parallel

#SLURM_ARRAY_TASK_ID=1

wd="/scratch/aob2x/daphnia_hwe_sims"
datadir="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase"
chr=$( grep "^${SLURM_ARRAY_TASK_ID}," ${datadir}/chrs.csv | cut -f2 -d',' )
echo $chr

### generate RABBIT input data
Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbitRedo/formatData.R ${chr}


### format RABBIT script file
echo "make mathematica input"
sed "s/STEM/${chr}/g" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbitRedo/template.m > \
${datadir}/${chr}.m

ls ${datadir}/${chr}.m


### Run RABBIT impute & reconstruct
math -script ${datadir}/${chr}.m

### convert paths
python /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/parseHaplotypes.py \
${chr}.all.csv \
${datadir}/ >> ${chr}.all.haps
