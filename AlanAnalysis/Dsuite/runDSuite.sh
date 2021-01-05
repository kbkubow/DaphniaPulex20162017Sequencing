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




### hybrid model
  module load bcftools
  bcftools view \
  -O v \
  -o /scratch/aob2x/daphnia_hwe_sims/dsuite/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.bcf


  sed -i 's/outgroup/Outgroup/g' /scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy.txt
  sed -i 's/Dcat/DCat/g' /scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy.txt
  sed -i 's/Dbunk/DBunk/g' /scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy.txt

  module load gcc/9.2.0

  cp /scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy.txt /scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy.3species.txt

  # nano /scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy.txt
  # cat /scratch/aob2x/daphnia_hwe_sims/dsuite/MapJune2020_ann.hyrbid_strategy.whatshap.shapeit.vcf | head -n1000 > /scratch/aob2x/daphnia_hwe_sims/dsuite/tmp.vcf
  /scratch/aob2x/daphnia_hwe_sims/dsuite/Dsuite/Build/Dsuite \
  Dtrios \
  /scratch/aob2x/daphnia_hwe_sims/dsuite/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf \
  /scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy.3species.txt

  cat hybrid_strategy.3species_BBAA.txt | awk '{if($3=="pulicaria") print $0}'






  nano /scratch/aob2x/daphnia_hwe_sims/dsuite/DBunk_D8_pulicaria.txt

  /scratch/aob2x/daphnia_hwe_sims/dsuite/Dsuite/Build/Dsuite \
  Dinvestigate \
  -w 100000,50000 \
  -n DBunk_D8_pulicaria \
  /scratch/aob2x/daphnia_hwe_sims/dsuite/MapJune2020_ann.hyrbid_strategy.whatshap.shapeit.vcf \
  /scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy.txt \
  /scratch/aob2x/daphnia_hwe_sims/dsuite/DBunk_D8_pulicaria.txt &
