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
