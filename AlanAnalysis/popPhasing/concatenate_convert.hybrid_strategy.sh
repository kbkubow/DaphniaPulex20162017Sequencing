module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822

ls -d /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/Scaffold*.whatshapp.onPerSC.hybrid_strategy.shapeit.bcf > /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/bcfs.hybrid_strategy.list

for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/Scaffold*.whatshapp.onPerSC.hybrid_strategy.shapeit.bcf; do
  #f=/scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/April17_2018_D8_Male1.Scaffold_1931_HRSCAF_2197.phase.vcf.gz
  echo "Index file -> $f"
  bcftools index \
  -f \
  ${f}
done

bcftools \
concat \
-f /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/bcfs.hybrid_strategy.list \
-O b \
-l \
-o /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.pulexOnly_hyrbidstrategy.whatshap.shapeit.bcf

bcftools \
view \
-O b \
-o /scratch/aob2x/daphnia_hwe_sims/popPhase/outgroup.bcf \
/scratch/aob2x/daphnia_hwe_sims/popPhase/outgroup.vcf

bcftools index /scratch/aob2x/daphnia_hwe_sims/popPhase/outgroup.bcf
bcftools index /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.pulexOnly_hyrbidstrategy.whatshap.shapeit.bcf

bcftools \
merge \
-O b \
-o /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.whatshap.shapeit.bcf \
/scratch/aob2x/daphnia_hwe_sims/popPhase/outgroup.bcf \
/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.pulexOnly_hyrbidstrategy.whatshap.shapeit.bcf



bcftools view \
/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.whatshap.shapeit.bcf |
grep -v "##" | less -S


bcftools view -h /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.whatshap.shapeit.bcf | \
grep -v "##" | cut -f10- | tr '\t' '\n'  | nl
