#!/usr/bin/env bash
#SBATCH -J ASE_readcounter_pul
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-4:00:00
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter_pul.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter_pul.%A_%a.err # Standard error

# ijob -c1 -p standard -A berglandlab
# submit as: sbatch --array=1-2 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/aseReadCounter_mf.sh


module load gatk/4.0.0.0

#cat /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf | awk '{
#  if(substr($1, 0, 1)=="#") {
#    print $0
#  } else {
#    for(i=1; i<=9; i++) printf $i"\t"
#    printf "0/1\t0/1\t0/1\t0/1\t0/1\t0/1\n"
#  }
#}' > /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.allvariant.vcf
#

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf
#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.allvariant.vcf

###

  gatk ASEReadCounter \
  --I /scratch/aob2x/daphnia_hwe_sims/pulicaria/pulicaria.sort.D84a.rmDup.bam \
  --variant /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.allvariant.vcf \
  --output /scratch/aob2x/daphnia_hwe_sims/pulicaria/pulicaria.sort.D84a.rmDup.aseReadCounter.allvariant.delim \
  --reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa
