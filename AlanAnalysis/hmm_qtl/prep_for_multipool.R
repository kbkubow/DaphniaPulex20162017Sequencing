### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(doMC)

### load data from Rivanna
  load("/nv/vol186/bergland-lab/alan/pool_AB.Rdata")
  m[!is.na(A),A.geno := unlist(sapply(m[!is.na(A)]$A, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[!is.na(B),B.geno := unlist(sapply(m[!is.na(B)]$B, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[,BA.geno.diff := B.geno - A.geno]
  m[,inform:=!((A.geno==0 & B.geno==0) | (A.geno==2 & B.geno==2))]


### generate weighted allele frequencies from pools
  m[,D8Male.f:=(effRD_D8Male1*effPA_D8Male1 + effRD_D8Male2*effPA_D8Male2) / (effRD_D8Male1+effRD_D8Male2)]
  m[,D8Male.rd:=(effRD_D8Male1+effRD_D8Male2)]

  m[,D8PE.f:=(effRD_D8PE1*effPA_D8PE1 + effRD_D8PE2*effPA_D8PE2) / (effRD_D8PE1+effRD_D8PE2)]
  m[,D8PE.rd:=(effRD_D8PE1+effRD_D8PE2)]

### polarize
  m[,A.geno.polarize := A.geno]
  m[,B.geno.polarize := B.geno]

  m[A.geno==2, A.geno.polarize := 0]
  m[A.geno==2 & B.geno==0, B.geno.polarize := 2]
  m[A.geno==2 & B.geno==2, B.geno.polarize := 0]
  m[,BA.geno.diff.polarize := B.geno.polarize - A.geno.polarize]

  m[,D8Male.f.polarize:=D8Male.f]
  m[,D8PE.f.polarize:=D8PE.f]
  m[A.geno==2, D8Male.f.polarize := 1 - D8Male.f]
  m[A.geno==2, D8PE.f.polarize := 1 - D8PE.f]

### export file
  m <- m[(D8Male.f.polarize/2 + D8PE.f.polarize/2)>.05 & (D8Male.f.polarize/2 + D8PE.f.polarize/2)<.95]
  m.ag <- m[,.N, chr]
  m <- m[chr%in%m.ag[N>1000]$chr]

  ex.male <- data.table(chr=m$chr, pos=m$pos, A=m$D8Male.f.polarize*m$D8Male.rd, B=(1-m$D8Male.f.polarize)*m$D8Male.rd)
  ex.pe <- data.table(chr=m$chr, pos=m$pos, A=m$D8PE.f.polarize*m$D8PE.rd, B=(1-m$D8PE.f.polarize)*m$D8PE.rd)

  foreach(i=unique(ex.male$chr))%do%{
    print(i)
    write.table(ex.male[chr==i, c("pos", "A", "B"),with=F], quote=F, row.names=F, col.names=F,
                file=paste("/scratch/aob2x/daphnia_hwe_sims/multipool/inputData/", i, ".male.count.delim", sep=""), sep="\t")

    write.table(ex.male[chr==i, c("pos", "A", "B"),with=F], quote=F, row.names=F, col.names=F,
                file=paste("/scratch/aob2x/daphnia_hwe_sims/multipool/inputData/", i, ".pe.count.delim", sep=""), sep="\t")

  }





#### with ASE read counter data
### load ASEReadCounter data
  D8Male <- fread("/nv/vol186/bergland-lab/alan/D8Male.pooledAF.aseReadCounter.allvariant.delim")
  D8PE <- fread("/nv/vol186/bergland-lab/alan/D8PE.pooledAF.aseReadCounter.allvariant.delim")

### merge
  setkey(D8Male, chr, pos)
  setkey(D8PE, chr, pos)

  m <- merge(D8Male[,c("chr", "pos", "altCount", "totalCount"), with=F], D8PE[,c("chr", "pos", "altCount", "totalCount"), with=F])
  m[,freq:=(altCount.x+altCount.y)/(totalCount.x+totalCount.y)]
  m <- m[freq>.05 & freq<.95]

  m.ag <- m[,.N,chr]
  m <- m[chr%in%m.ag[N>1000]$chr]
  m[,refCount.x:=totalCount.x - altCount.x]
  m[,refCount.y:=totalCount.y - altCount.y]

### output
  foreach(i=unique(m$chr))%do%{
    #i<-unique(m$chr)[1]
    print(i)
    write.table(m[chr==i, c("pos", "altCount.x", "refCount.x"),with=F], quote=F, row.names=F, col.names=F,
                file=paste("/scratch/aob2x/daphnia_hwe_sims/multipool/inputData/", i, ".male.count.delim", sep=""), sep="\t")

    write.table(m[chr==i, c("pos", "altCount.y", "refCount.y"),with=F], quote=F, row.names=F, col.names=F,
                file=paste("/scratch/aob2x/daphnia_hwe_sims/multipool/inputData/", i, ".pe.count.delim", sep=""), sep="\t")

  }















### load data from Rivanna
  load("/nv/vol186/bergland-lab/alan/pool_AB.Rdata")
  m[!is.na(A),A.geno := unlist(sapply(m[!is.na(A)]$A, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[!is.na(B),B.geno := unlist(sapply(m[!is.na(B)]$B, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[,BA.geno.diff := B.geno - A.geno]
  m[,inform:=!((A.geno==0 & B.geno==0) | (A.geno==2 & B.geno==2))]


### merge
  setnames(D8Male, c("contig", "position"), c("chr", "pos"))
  setnames(D8PE, c("contig", "position"), c("chr", "pos"))


  ### male pool
    setkey(D8Male, chr, pos)
    D8Male[,altFreq:=altCount/totalCount]




    male <- D8Male[altFreq>.05 & altFreq<.95]
    male[,refCount:=totalCount - altCount]

  ### PE pool
    setkey(D8PE, chr, pos)
    D8PE[,altFreq:=altCount/totalCount]
    pe <- D8PE[altFreq>.05 & altFreq<.95]
    pe[,refCount:=totalCount - altCount]

  ### output
    foreach(i=unique(ex.male$chr))%do%{
      print(i)
      write.table(ex.male[chr==i, c("pos", "A", "B"),with=F], quote=F, row.names=F, col.names=F,
                  file=paste("/scratch/aob2x/daphnia_hwe_sims/multipool/inputData/", i, ".male.count.delim", sep=""), sep="\t")

      write.table(ex.male[chr==i, c("pos", "A", "B"),with=F], quote=F, row.names=F, col.names=F,
                  file=paste("/scratch/aob2x/daphnia_hwe_sims/multipool/inputData/", i, ".pe.count.delim", sep=""), sep="\t")

    }
