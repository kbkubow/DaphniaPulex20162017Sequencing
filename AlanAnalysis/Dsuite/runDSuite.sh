module load gcc/9.2.0

/scratch/aob2x/daphnia_hwe_sims/dsuite/Dsuite/Build/Dsuite \
Dtrios \
/scratch/aob2x/daphnia_hwe_sims/dsuite/daphnids.orig.vcf \
/scratch/aob2x/daphnia_hwe_sims/dsuite/SETS.txt



/scratch/aob2x/daphnia_hwe_sims/dsuite/Dsuite/Build/Dsuite \
Dinvestigate \
-w 100000,50000 \
-n D8_Dbunk_pulicaria \
/scratch/aob2x/daphnia_hwe_sims/dsuite/daphnids.orig.vcf \
/scratch/aob2x/daphnia_hwe_sims/dsuite/SETS.txt \
/scratch/aob2x/daphnia_hwe_sims/dsuite/test_trios.txt

module load bcftools
bcftools view 
module load gcc/9.2.0

/scratch/aob2x/daphnia_hwe_sims/dsuite/Dsuite/Build/Dsuite \
Dtrios \
/scratch/aob2x/daphnia_hwe_sims/dsuite/daphnids.orig.vcf \
/scratch/aob2x/daphnia_hwe_sims/dsuite/SETS.txt
