#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem 24G
#SBATCH -t 0-20:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

# submit as
# sbatch --array=1-$( wc -l /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f1 -d' ' ) doHarp.slurm

module load gparallel
module load gcc/7.1.0  openmpi/3.1.4 boost/1.60.0

# copy over bam files to scratch and rename index
  # cp /project/berglandlab/Karen/MappingDec2019/MaleFemale2018Pools/bams/*D8* \
  # /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/.

  # mv /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8Male1.filt.merged.mdup.bai \
  # /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8Male1.filt.merged.mdup.bam.bai \

  # mv /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8Male2.filt.merged.mdup.bai \
  # /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8Male2.filt.merged.mdup.bam.bai \

  # mv /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8PE1.filt.merged.mdup.bai \
  # /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8PE1.filt.merged.mdup.bam.bai \

  # mv /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8PE2.filt.merged.mdup.bai \
  # /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8PE2.filt.merged.mdup.bam.bai \


# get chromosome from job id
  #SLURM_ARRAY_TASK_ID=1
  chromosome=$( grep "^"${SLURM_ARRAY_TASK_ID}" " /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' )
  chrLength=$( grep "^"${SLURM_ARRAY_TASK_ID}" " /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f3 -d' ' )
  bam=$( grep "^"${SLURM_ARRAY_TASK_ID}" " /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f4 -d' ' )
  sample=$( echo $bam | rev  | cut -f1 -d'/' | rev | cut -f1 -d'.' | cut -f3 -d"_" )

# harp parameters
  window_step=100000
  window_width=100000
  #window_step=100
  #window_width=100

# paths
  harp="/project/berglandlab/cory/2019-10-25/genome-reconstruction-revision/04_RABBIT/harp"
  referenceGenome="/project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa"
  bamFolder="/scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/"
  priorsFile="/scratch/aob2x/daphnia_hwe_sims/harp_pools/priors/${chromosome}.csv"

# Run harp like
  echo running harp like
  $harp like \
  --bam ${bam} \
  --region ${chromosome}:1-${chrLength} \
  --refseq ${referenceGenome} \
  --snps ${priorsFile} \
  --stem /scratch/aob2x/daphnia_hwe_sims/harp_pools/out/${sample}/${sample} \
  -v \
  --out /scratch/aob2x/daphnia_hwe_sims/harp_pools/out/${sample}

  echo running harp freq
  $harp freq \
  --bam ${bam} \
  --region $chromosome:1-${chrLength} \
  --refseq $referenceGenome \
  --snps ${priorsFile} \
  --stem /scratch/aob2x/daphnia_hwe_sims/harp_pools/out/${sample}/${sample} \
  --out /scratch/aob2x/daphnia_hwe_sims/harp_pools/out/${sample}/  \
  --window_step $window_step \
  --window_width $window_width \
  --em_min_freq_cutoff 0.0001
