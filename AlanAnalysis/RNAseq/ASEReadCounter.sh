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

# submit as: sbatch --array=1-8 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/ASEReadCounter.sh
# sacct -j 20416102
# cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.17863582_1.out

module load gatk/4.0.0.0

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf
#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.allvariant.vcf

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf

#SLURM_ARRAY_TASK_ID=2
wd=/scratch/aob2x/daphnia_hwe_sims/

samp=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/samples )
echo $samp


###

gatk ASEReadCounter \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}.trim.bam \
--variant /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf \
--output /scratch/aob2x/daphnia_hwe_sims/rnaseq/ase/${samp}_rnaseq_asereadcounter.delim \
--reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa
