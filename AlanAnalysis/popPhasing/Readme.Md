1. DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_prepareSplitVCF.sh
2. DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_SplitVCF.sh



A few notes

### update bai files (I think they must get copied before the bams which throws the error re: their age)
  #for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/bams/*.bai; do
  #  touch -c ${f}
  #done

### one of the bams is problematic, remove from analysis
  #samtools quickcheck -v /scratch/aob2x/daphnia_hwe_sims/popPhase/bams/*.bam > ~/bad_bams.fofn
  #
  #samtools quickcheck /project/berglandlab/Karen/MappingDec2019/bams/PulexBams/Spring_2017_DBunk_340_finalmap_mdup.bam
  #
  #grep -m1 "#CHROM" /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf | \
  #tr '\t' '\n' | grep "Spring_2017_DBunk_340"

### Karen says that these are the same individual
  #Spring_2016_D8_8.1
  #Spring_2016_D8_8.10

### get 12 big chromosomes
#cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' \
#| sort | uniq | grep -v "chr" | sed 's/^/--chromosome /g' | sed 's/$/ \\/g'
