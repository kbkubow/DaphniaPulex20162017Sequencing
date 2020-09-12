#!/usr/bin/env bash
#
#SBATCH -J trioPhase_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 1:00:00 # Running time of 1 hours
#SBATCH --mem 24G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run as
# sbatch --array=1-12 ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.runRabbit.sh
# sacct -j 15041987
# cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.14761054_3.out

# sbatch --array=1 ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.runRabbit.sh
# sacct -j 14874181
# cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.14874181_1.out



module load intel/18.0 intelmpi/18.0 R/3.6.3
module load mathematica
module load parallel

#SLURM_ARRAY_TASK_ID=11
cm=10
wd="/scratch/aob2x/daphnia_hwe_sims"
datadir=/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm
chr=$( grep "^${SLURM_ARRAY_TASK_ID}," ${datadir}/chrs.csv | cut -f2 -d',' )
echo $chr $datadir

### make dir
mkdir -p ${wd}/Rabbit_phase_10cm/${chr}

### generate RABBIT input data
### which f1s
  # options: onlyPheno_AxC; wildF1s_AxC; all_AxC; all_CxC
  set="all_AxC"

    Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.formatData.R ${chr} ${cm} ${set}


### format RABBIT script file
echo "make mathematica input"
sed "s/STEM/${chr}/g" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/template.m > \
${datadir}/${chr}/${chr}.m

ls ${datadir}/${chr}/${chr}.m


### Run RABBIT impute & reconstruct
math -script ${datadir}/${chr}/${chr}.m

### convert paths
# datadir=/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm_old
python ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/old/rabbit.parseHaplotypes.py \
${chr}.all.csv \
${datadir}/${chr}/ >> ${datadir}/${chr}/${chr}.haps.delim
