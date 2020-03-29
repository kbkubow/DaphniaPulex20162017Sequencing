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

### set filter
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=snps.use$id, sample.id=sc[Species=="pulex"]$clone)

### get King Relatedness estimate
  king <- snpgdsIBDKING(genofile, sample.id=sc[Species=="pulex"]$clone, snp.id=snps.use$id, autosome.only=FALSE,
                remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
                type=c("KING-robust"),
                family.id=NULL, num.thread=10L, useMatrix=FALSE, verbose=TRUE)

  king.rel <- king$kinship
  dimnames(king.rel) <- list(king$sample.id, king$sample.id)

  king.dt <- as.data.table(melt(king.rel))
  setnames(king.dt, c("Var1", "Var2"), c("ID1", "ID2"))

### import truffle data
  truffle <- fread("/scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd.ibd")

### merge
  setkey(king.dt, ID1, ID2)
  setkey(truffle, ID1, ID2)

  m <- merge(king.dt, truffle)

### format
  " the input file must have columns 1-4 and 7-10 described in option --plink_ibd, except that the PI_HAT/RELATEDNESS column can be a different measure of relatednes"
  "FID1(1) IID1(2) FID2(3) IID2(4) RT(5) EZ(6) IDB0(7) IBD1(8) IBD2(9) PI_HAT(10)"

  out <- data.table(FID1=m$ID1, IID1=m$ID1, FID2=m$ID2, IID2=m$ID2, RT=0, EZ=0, IBD0=m$IBD0, IBD1=m$IBD1, IBD2=m$IBD2, PI_HAT=m$value)

  write.table(out, file="/scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim",
              quote=F, row.names=F, sep="\t")
