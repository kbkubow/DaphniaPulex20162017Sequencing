#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-02:00:00 # Running time of 4 days
#SBATCH --mem 100G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/maketree.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/maketree.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# ijob -c1 -p standard -A berglandlab
### run with: sbatch --array=1-12 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/combineFASTA.daphnids.sh
### sacct -u aob2x -j 19110776
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.19110025_10.err

# sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.makeTrees.sh
# sacct -j 19141305

### modules
  module load samtools parallel

### define parameters
  chr=Scaffold_2217_HRSCAF_2652
  start=5173221
  stop=5223221


## set up RAM disk
  ## rm /scratch/aob2x/test/*
  #tmpdir="/scratch/aob2x/test"
  #SLURM_JOB_ID=1; SLURM_ARRAY_TASK_ID=1

  [ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

  tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

### make bed file
  echo ${chr}","${start}","${stop} | tr ',' '\t' > ${tmpdir}/region.bed

### extract chromosome from each individual
  getRegion () {
    chr=${1}
    file=${2} # file="/scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/Spring_2017_DBunk_360.1.fa"
    tmpdir=${3}

    stem=$( echo ${file} | rev | cut -d'/' -f1 | rev )

    echo ${chr} ${stem}

    if [ ! -f ${file}.fai ]; then samtools faidx ${file}; fi

    ~/seqtk/seqtk \
    subseq \
    ${file} \
    ${tmpdir}/region.bed | sed "s/${chr}/${stem};${chr}/g" > ${tmpdir}/${chr}_${stem}

  }
  export -f getRegion

  parallel getRegion ::: ${chr} ::: $( ls -d  /scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/*.fa ) ::: ${tmpdir}

### combine
  awk 'NR>1 && FNR==1{print ""};1' ${tmpdir}/*.fa > ${tmpdir}/${chr}_${start}_${stop}.fasta

### make tree
  /home/aob2x/iqtree-1.6.12-Linux/bin/iqtree \
  -nt 20 \
  -redo \
  -s ${tmpdir}/${chr}_${start}_${stop}.fasta

  cp ${tmpdir}/${chr}_${start}_${stop}.fasta* /scratch/aob2x/daphnia_hwe_sims/popPhase/trees/

  rm -fr ${tmpdir}

##
