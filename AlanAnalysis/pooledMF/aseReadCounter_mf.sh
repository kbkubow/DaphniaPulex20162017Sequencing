#!/usr/bin/env bash
#SBATCH -J ASE_readcounter
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-4:00:00
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.%A_%a.err # Standard error

# ijob -c1 -p standard -A berglandlab
# submit as: sbatch --array=1-2 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/aseReadCounter_mf.sh


module load gatk/4.0.0.0

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf

###

if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]; then

  gatk ASEReadCounter \
  --I /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8PE1.filt.merged.mdup.bam \
  --I /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8PE2.filt.merged.mdup.bam \
  --variant /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf \
  --output /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8PE.pooledAF.aseReadCounter.delim \
  --reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa \

elif [[ ${SLURM_ARRAY_TASK_ID} -eq 2 ]]; then

  gatk ASEReadCounter \
  --I /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8Male1.filt.merged.mdup.bam \
  --I /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8Male2.filt.merged.mdup.bam \
  --variant /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf \
  --output /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8Male.pooledAF.aseReadCounter.delim \
  --reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa \

fi


SLURM_ARRAY_TASK_ID=1
if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]
then
  echo "foo"
elif [[ ${SLURM_ARRAY_TASK_ID} -eq 2 ]]
then
  echo "bar"
fi
