#!/usr/bin/env bash
#SBATCH -J ASE_readcounter
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-1:00:00
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.%A_%a.err # Standard error

# ijob -c1 -p standard -A berglandlab
# submit as: sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/aseReadCounter_mf.sh


module load gatk/4.0.0.0

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.vcf

gatk ASEReadCounter \
--I /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8PE1.filt.merged.mdup.bam \
--I /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8PE2.filt.merged.mdup.bam \
--I /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8Male1.filt.merged.mdup.bam \
--I /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/HT2LNDSXX_s1_D8Male2.filt.merged.mdup.bam \
--variant /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.vcf \
--output /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/pooledAF.aseReadCounter.delim \
--reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa


mkdir /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/
