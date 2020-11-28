# A. Convert original VCF to BCF and filter on sites that pass quality filtering via `snpsvarpulexpresentinhalf_table_20200623`
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/filterVCF.sh`

# B. Make jobs file (1 line per chr X clone, downsampled to 1 clone/SC, no lab F1s)
  1. `onePerSC.R` this generates 1 file which is the 1 line per chr per SC (passing filter)

# C. Run Whatshapp to generate read-backed phase-sets
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_SplitVCF.sh`

# D. merge individual vcfs
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/mergeVCF.sh`

# E. shapeit4
  ### Run shapeit on each chromosome separatley
    1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/shapeit4.sh`

  ### concatenate resulting VCFs
    2.

# F. KING to get basic relationships; is based on observed (not imputed) data; identifies parent offspring pairs
  ### concatenate whatshapp VCF files, rename chromosomes, convert bcf to plink
    1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/bcf2plink.sh`

  ### download KING
    cd /scratch/aob2x/daphnia_hwe_sims/king
    wget http://people.virginia.edu/~wc9c/KING/Linux-king.tar.gz
    tar xzvf Linux-king.tar.gz

  ### run KING
    2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/king.sh`

# G. duoHMM
  ## fixes phasing errors using pedigree; generates recombination map file


  ### poking about
    module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822
    bcftools view Scaffold_2373_HRSCAF_2879.whatshapp.shapeit.bcf | grep -v "##" | less -S

# duoHMM
  wget https://mathgen.stats.ox.ac.uk/genetics_software/duohmm/duohmm_v0.1.7.tar.gz
  wget https://github.com/molgenis/systemsgenetics/releases/download/1.4.0_20-8.1/GenotypeHarmonizer-1.4.23-dist.tar.gz

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
