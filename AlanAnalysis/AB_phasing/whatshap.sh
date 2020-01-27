#!/usr/bin/env bash
#
#
#SBATCH -J trioPhase_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=8 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:90 # Running time of 90 minutes
#SBATCH --mem 8G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#ijob -c1 -p standard -A berglandlab
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools

### install whatsapp
  #conda create -y -n whatshapp-env

  #source activate whatshapp-env
  #conda config --add channels bioconda

  #conda install whatshap nomkl

### load env
  source activate whatshapp-env
  whatshap_dir=/home/aob2x/.conda/pkgs/whatshap-0.18-py36h6bb024c_0/bin

### run whatshap

  /home/aob2x/.conda/envs/whatshapp-env/bin/whatshap \
  phase \
  --ped /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.ped \
  -o /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.phased.vcf \
  /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.vcf

  #\
  #/scratch/aob2x/daphnia_hwe_sims/trioPhase/April17_2018_D8_Male1_finalmap_mdup.bam \
  #/scratch/aob2x/daphnia_hwe_sims/trioPhase/April_2017_D8_103_finalmap_mdup.bam \
  #/scratch/aob2x/daphnia_hwe_sims/trioPhase/bams/April_2017_D8_125_finalmap_mdup.bam


  File "/home/aob2x/.conda/envs/whatshapp-env/lib/python3.6/site-packages/whatshap/utils.py", line 37, in detect_file_format
