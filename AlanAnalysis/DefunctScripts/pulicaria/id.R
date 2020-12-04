### ijob --mem=24G -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R


### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### load North American pulicaria
  puli <- fread("/scratch/aob2x/daphnia_hwe_sims/pulicaria/pulicaria.sort.D84a.rmDup.aseReadCounter.allvariant.delim")
  setnames(puli, c("contig", "position"), c("chr", "pos"))
  setkey(puli, chr, pos)
  puli[,puli.geno:=round(refCount/totalCount*10)/10]
  puli[puli.geno!=1 & puli.geno!=0, puli.geno:=.5]

### open GDS file
  genofile <- seqOpen("/scratch/aob2x/daphnia_hwe_sims/pulicaria/Males_2018_filtsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"), pos=seqGetData(genofile, "position"), id=seqGetData(genofile, "variant.id"))

### merge
  setkey(puli, chr, pos)
  setkey(snp.dt, chr, pos)

  snp.dt <- merge(snp.dt, puli[puli.geno%in%c(0,0.5,1),c("chr", "pos", "puli.geno"), with=F])
  snp.dt[,puli.geno:=2*puli.geno]

### get genotypes
  seqSetFilter(genofile, variant.id=snp.dt$id)
  genomat <- seqGetData(genofile, "$dosage")

### compare
  o <- foreach(i=1:dim(genomat)[1], .combine="rbind")%do%{
    print(i)
    data.table(i=i, TT=sum(genomat[i,]==snp.dt$puli.geno, na.rm=T), FF=sum(genomat[i,]!=snp.dt$puli.geno, na.rm=T))
  }

  o[,frac:=TT/(TT+FF)]
  o[,N:=TT+FF]

  pulicaria <- seqGetData(genofile, "sample.id")[o[N>1000][frac>.9]$i]

### save british pulicaria genome data
  seqResetFilter(genofile)
  seqSetFilter(genofile, sample.id=pulicaria)

  pu.dt <- as.data.table(t(seqGetData(genofile, "$dosage")))
  pu.dt[,chr:=seqGetData(genofile, "chromosome")]
  pu.dt[,pos:=seqGetData(genofile, "position")]

  setnames(pu.dt, paste("V", c(1:5), sep=""), pulicaria)

  pu.dt.l <- melt(pu.dt, id.vars=c("chr", "pos"))

  pu.ag <- pu.dt.l[,list(mean.geno=mean(value, na.rm=T), n=sum(!is.na(value))), list(chr, pos)]
  pu.ag[,puli.geno:=mean.geno]
  pu.ag[!mean.geno%in%c(0,1,2) & n>0, puli.geno:=1]

  prop.table(table(pu.ag[n>0][puli.geno%in%c(0,1,2)]$puli.geno))

### save
  save(pu.ag, file="/scratch/aob2x/daphnia_hwe_sims/pulicaria/puAg.Rdata")
  save(pu.ag, file="/nv/vol186/bergland-lab/alan/puAg.Rdata")
