#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

### libraries
  library(SeqArray)
  library(data.table)
  #library(SNPRelate)


### open GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")

### clone metadata file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/CloneInfoFilePulexandObtusa_withmedrd_update20200324")
  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]

  sc.ag <- sc[Species=="pulex",list(clone=clone[which.max(medrd)][1]), list(SC.uniq)]
  table(sc.ag$SC.uniq)

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

### export
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=snps.use$id, sample.id=sc.ag$clone)

  seqGDS2VCF(genofile,
             "/scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.vcf.gz",
             info.var=character(), fmt.var=character(), use_Rsamtools=TRUE, verbose=TRUE)


  zcat /scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.vcf.gz | grep -v "##" | head -n10 | awk '{print NF}'

### rename chromosomes
  rename.chr <- snps.ag[N>5000]
  rename.chr[,newChr:=c(1:12)]

  write.table(rename.chr[,c("chr", "newChr"),with=F], file="/scratch/aob2x/daphnia_hwe_sims/pedigree/chr_rename.delim",
              quote=F, row.names=F, col.names=F)
