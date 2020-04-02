#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

### libraries
  library(data.table)
  #library(GWASTools)
  library(SeqArray)
  library(SNPRelate)
  source("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/ibdAssignKing.R")
  
### open GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")

### clone metadata file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/CloneInfoFilePulexandObtusa_withmedrd_update20200324")

  sc <- sc[Nonindependent==0]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]

### hard filtering of SC
  #sc <- sc[SC!="OO"]

  sc.ag <- sc[Species=="pulex", list(clone=clone[which.max(medrd)][1], year=year[which.min(year)][1]), list(SC.uniq=paste(SC.uniq, clone, sep=";"))]
  sc.ag[,age:=max(sc.ag$year) - year + 1]
  table(sc.ag$SC.uniq)


### filtered SNP set, LD-pruned by Karen
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
    seqSetFilter(genofile, variant.id=snps.use$id, sample.id=sc.ag$clone)

  ### get King Relatedness estimate
    king <- snpgdsIBDKING(genofile, snp.id=snps.use$id, autosome.only=FALSE,
                  remove.monosnp=TRUE, maf=.01, missing.rate=NaN,
                  type=c("KING-robust"),
                  family.id=NULL, num.thread=10L, useMatrix=FALSE, verbose=TRUE)

                  king <- snpgdsIBDSelection(king, kinship.cutoff=1/32)


### relatedness




o <- ibdAssignRelatednessKing(ibs0=king$IBS0, kc=king$kinship, cut.kc.dup=1/(2^(3/2)),
                       cut.kc.fs=1/(2^(5/2)), cut.kc.deg2=1/(2^(7/2)),
                       cut.kc.deg3=1/(2^(9/2)), cut.ibs0.err=0.003)

ibdPlot(k0=ibd$IBD0, k1=ibd$IBD1, relation=ibdAssignRelatedness(k0=ibd$IBD0, k1=ibd$IBD1))


, relation=NULL, color=NULL,
        rel.lwd=2, rel.draw=c("FS", "Deg2", "Deg3"), ...)
