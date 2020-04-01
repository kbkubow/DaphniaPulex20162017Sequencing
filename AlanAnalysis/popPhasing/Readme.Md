# Get data in right place
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/copyBams.sh`

# these are some index files created by scripts for slurm. I don't remember exactly which scripts make these files.
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/bams2phase.sh`
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/path2bams`

# Run Whatshapp to generate read-backed phase-sets
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_SplitVCF.sh`

# merge individual vcfs
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/mergeVCF.sh`

# shapeit4
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/shapeit4.sh`

  ### poking about
    module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822
    bcftools view Scaffold_2373_HRSCAF_2879.whatshapp.shapeit.bcf | grep -v "##" | less -S

# duoHMM
  wget https://mathgen.stats.ox.ac.uk/genetics_software/duohmm/duohmm_v0.1.7.tar.gz

  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/duoHMM.sh`






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
