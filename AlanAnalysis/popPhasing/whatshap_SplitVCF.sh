#!/usr/bin/env bash
#
#
#SBATCH -J popPhasing_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0:15:00 # Running time of 15 minutes
#SBATCH --mem 4G # Memory request of 4GGB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# ijob -c1 -p standard -A berglandlab
### nJobs=$( wc -l /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.delim | cut -f1 -d' ' )
### run with: sbatch --array=1-${nJobs} /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_SplitVCF.sh
### sacct -u aob2x -j 10007470
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_whatshapp.10007470_20.out



### load modules
  module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822
  export PATH=$HOME/.local/bin:$PATH

### make job file
  #chr=$( cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' | sort | uniq | grep -v "chr" )

  makeJobs () {
    #line="foo"
    line=${1}

    grep -m1 "#CHROM" /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf |
    cut -f 10- | \
    tr '\t' '\n' | \
    sed "s/^/${line}\t/g"

  }
  export -f makeJobs

  #parallel -j 1 makeJobs ::: ${chr} > /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.delim
  #cat /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.delim | awk '{print NR"\t"$0}' > /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.delim


### get parameters
  # SLURM_ARRAY_TASK_ID=2
  chr=$( awk -v samp=${SLURM_ARRAY_TASK_ID} '{if(NR==samp) { print $2} }' < /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.delim )
  sample=$( awk -v samp=${SLURM_ARRAY_TASK_ID} '{if(NR==samp) { print $3} }' < /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.delim  )


  echo "chr: "${chr}
  echo "samp: "${sample}

### extract simplified vcf file per chromosome per sample
  echo "Getting simple vcf"

  bcftools view \
  -s ${sample} \
  -O v \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.${chr}.vcf.gz | \
  awk '{
    a=0
    if(substr($0, 0, 1)=="#") {
      print $0
    } else {
      printf $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t.\tGT\t"
      split($10, sp, ":")
      printf sp[1]"\n"
    }
  }' > /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${sample}_${chr}.vcf

### extract reads from bam file per chromosome per sample
  echo "getting bam"

  samtools view \
  -b \
  -q 20 \
  -M \
  -O BAM \
  -@ 1 \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/bams/${sample}_finalmap_mdup.bam \
  ${chr} > \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${sample}_${chr}.bam

  samtools index /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${sample}_${chr}.bam

### run whatshap
  echo "whatshapp"

  whatshap \
  phase \
  -o /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${sample}.${chr}.phase.vcf \
  --chromosome ${chr} \
  --sample ${sample} \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${sample}_${chr}.vcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${sample}_${chr}.bam

### clean up
  rm /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${sample}_${chr}.vcf
  rm /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${sample}_${chr}.bam
