### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(ggtern)
  library(SeqArray)
  library(doMC)
  registerDoMC(20)
  library(SeqVarTools)

### load precomputed file (HWE_simulations.prep.rivanna.R makes this file)
  load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
  sc[,year:=tstrsplit(clone, "_")[[2]]]
  sc[,pond:=gsub("Dcat", "DCat", pond)]

### open genotype file
  genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### get consensus sequences for list of clones
  consSeq <- function(clones, name) {
    #clones <- sc[SC=="A"]$clone; name="A"
    seqSetFilter(genofile, variant.id=snp.dt[(final.use)]$id, sample.id=clones)

    tmp <- data.table(af=seqAlleleFreq(genofile))
    tmp[,af.r:=round(af, 1)]
    tmp <- cbind(snp.dt[(final.use)], tmp)
    tmp[,name:=name]
  }

  A.cons <- consSeq(clones=sc[SC=="A"]$clone, name="A")
  B.cons <- consSeq(clones=sc[SC=="B"]$clone, name="B")

### merge
  setkey(A.cons, id)
  setkey(B.cons, id)

  m <- merge(A.cons[,c("id", "af.r"), with=T], B.cons[,c("id", "af.r"), with=T])

### ID sites
  #sites <- m[af.r.x==.5 & af.r.y==1]
  sites <- m[af.r.x%in%c(0, .5, 1) & af.r.y%in%c(0, .5, 1)]

### get stats
  seqResetFilter(genofile)
  #seqSetFilter(genofile, variant.id=sites$id, sample.id=sc[SC!="A" & SC!="B" & pond%in%c("D8", "DBunk")]$clone)
  seqSetFilter(genofile, variant.id=sites$id, sample.id=sc[SC!="A" & SC!="B" & pond%in%c("D8", "DBunk")]$clone)

  o <- as.data.table(hwe(genofile))

### do test
  sites[,zw_AB := af.r.x==.5 & af.r.y==1]

  o[,zw_pond:=NA]
  o[afreq<.95, zw_pond := naa<=1]

  m <- merge(sites, o, by.x="id", by.y="variant.id")

  fisher.test(table(m$zw_AB, m$zw_pond))

  table(m$zw_AB, m$zw_pond)

  boxplot(afreq~interaction(zw_AB, zw_pond), m)
  ggplot(data=m[zw_pond==T], aes(f, group=as.factor(zw_AB), color=as.factor(zw_AB))) + geom_density()


### sliding window
  setkey(snp.dt, id)
  setkey(m, id)

  mm <- merge(snp.dt, m)[!is.na(zw_AB) & !is.na(zw_pond)]
  mm[,posBin:=floor(pos/1e5)*1e5]

  mm.ag <- mm[,list(frac_zw_AB=mean(zw_AB), frac_zw_pond=mean(zw_pond), frac_both=mean(zw_AB & zw_pond), n=length(zw_AB)),
               list(chr, posBin)]

  mm.ag[,exp:=frac_zw_pond * frac_zw_AB]
  mm.ag[,en:=log2(frac_both/exp)]

  ggplot(data=mm.ag[n>50], aes(y=frac_both, x=posBin)) + geom_point() + facet_wrap(~chr)
  
  ggplot(data=mm.ag[n>50], aes(x=frac_zw_AB, y=frac_zw_pond, color=chr)) + geom_point()

  hist(m[zw_AB & zw_pond]$id, breaks=1000)


mm.ag[n>50][frac_both>.5]
