#!/usr/bin/env bash
#
#
#SBATCH -J popPhasing_mergeVCF # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 04:00:00 # Running time of 4 days
#SBATCH --mem 100G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run primus

  ### download
  wget https://primus.gs.washington.edu/docroot/versions/PRIMUS_v1.9.0.tgz


  ###
  wd=/scratch/aob2x/daphnia_hwe_sims/pedigree/PRIMUS_v1.9.0/bin


  module load perl

  ${wd}/run_PRIMUS.pl \
  -i FILE=/scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim \
  IBD0=7 IBD1=8 IBD2=9 PI_HAT=10 \
  -o /scratch/aob2x/daphnia_hwe_sims/pedigree/primusOut/
  -v 3 &



  scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/pedigree/primusOut/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim_network1.dot \
  .


    scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/pedigree/primusOut/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim_network2.dot \
    .
