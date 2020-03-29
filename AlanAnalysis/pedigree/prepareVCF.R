#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

### libraries
  library(SeqArray)
  library(data.table)
  library(SNPRelate)


### open GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")

### clone metadata file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/CloneInfoFilePulexandObtusa_withmedrd_update20200324")

### filtered SNP set
  snps <- fread("/project/berglandlab/Karen/MappingDec2019/finalsetsnpset01pulex_table_20200207")
  snps.ag <- snps[,list(.N), chr]

  snps[,goodChr:=F]
  snps[chr%in%snps.ag[N>5000]$chr, goodChr:=T]

### check
  seqResetFilter(genofile)
  seqSetFilterPos(genofile, chr=snps[goodChr==T]$chr, pos=snps[goodChr==T]$pos)

  snps.use <- data.table(chr=seqGetData(genofile, "chromosome"),
                           pos=seqGetData(genofile, "position"),
                          id=seqGetData(genofile, "variant.id"))

### LD-prune
  snpset <- snpgdsLDpruning(genofile,
                  sample.id = sc[Species=="pulex"]$clone,
                  snp.id = snps.use$id, autosome.only = FALSE, slide.max.bp = 5000,
                  slide.max.n = NA, ld.threshold = 0.2, num.thread = 1, verbose = TRUE)

  snpset <- unlist(snpset)


### export
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=snpset, sample.id=sc[Species=="pulex"]$clone)

  seqGDS2VCF(genofile, 
             "/scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.vcf.gz",
             info.var=character(), fmt.var=character(), use_Rsamtools=TRUE, verbose=TRUE)
