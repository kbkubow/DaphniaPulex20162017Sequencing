#!/usr/bin/env bash
#
#
#SBATCH -J trioPhase_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 12:00:00 # Running time of 12 hours
#SBATCH --mem 12G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#ijob -c1 -p standard -A berglandlab
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools
export PATH=$HOME/.local/bin:$PATH

### run whatshap

whatshap \
phase \
--ped /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.ped \
-o /scratch/aob2x/daphnia_hwe_sims/trioPhase/MapDec19PulexOnlyB_filtsnps10bpindels_snps_filter_pass_lowGQmiss.trioPhased.vcf \
/project/berglandlab/Karen/MappingDec2019/MapDec19PulexOnlyB_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf
