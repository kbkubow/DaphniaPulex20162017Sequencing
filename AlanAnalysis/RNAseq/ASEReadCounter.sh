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

# submit as: sbatch --array=1-4 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/aseReadCounter_mf.sh
# sacct -j 17864732
# cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.17863582_1.out

module load gatk/4.0.0.0

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf
#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.allvariant.vcf

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf


###

gatk ASEReadCounter \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/d8_179_1.trim.bam \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/d8_179_2.trim.bam \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/d8_222_1.trim.bam \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/d8_222_2.trim.bam \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/d8_349_1.trim.bam \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/d8_349_2.trim.bam \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/d8_515_1.trim.bam \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/d8_515_2.trim.bam \
--variant /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf \
--output /scratch/aob2x/daphnia_hwe_sims/rnaseq/rnaseq_asereadcounter.delim \
--reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa
