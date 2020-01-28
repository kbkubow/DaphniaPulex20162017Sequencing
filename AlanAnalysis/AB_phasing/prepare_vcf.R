#ijob -c1 -p standard -A berglandlab
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### load SuperClone & SNP filter file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/Superclones201617182019pulexonlyD82016problematic_20200122")
  load("/project/berglandlab/Karen/MappingDec2019/snpsvarpulexpresentinhalf_20200121.Rdata")

### open GDS
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/MapDec19PulexOnlyB_filtsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"), pos=seqGetData(genofile, "position"), id=seqGetData(genofile, "variant.id"), key="id")
  use <- data.table(id=snpsvarpulexpresentinhalf, pass=T, key="id")
  snp.dt <- merge(snp.dt, use)

### 1. identify fixed difference between A&B
  ab.fd <- foreach(sc.i=c("A", "B"), .combine="cbind")%do%{
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sc[SC==sc.i]$clone, variant.id=snp.dt$id)

    data.table(af=seqAlleleFreq(genofile))
  }
  setnames(ab.fd, c(1,2), c("af.A", "af.B"))
  ab.fd <- cbind(ab.fd, snp.dt)
  ab.fd[!is.na(af.A),A.geno := unlist(sapply(ab.fd[!is.na(af.A)]$af.A, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  ab.fd[!is.na(af.B),B.geno := unlist(sapply(ab.fd[!is.na(af.B)]$af.B, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]


### 2. Identify F1s
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=snp.dt$id)

  genomat <- seqGetData(genofile, "$dosage")

  het.count <- foreach(i=1:dim(genomat)[1], .combine="rbind")%do%{
    print(paste(i, dim(genomat)[1], sep=" / "))

    tmp <- ab.fd
    tmp[,geno:=genomat[i,]]
    nLoci <- dim(tmp[!is.na(A.geno) & !is.na(B.geno) & !is.na(geno)])[1]
    tmp.ag <- tmp[!is.na(A.geno) & !is.na(B.geno),list(nRR=sum(geno==2, na.rm=T),
                                                      nRA=sum(geno==1, na.rm=T),
                                                      nAA=sum(geno==0, na.rm=T),
                                                      n=sum(!is.na(geno))),
                                                  list(A.geno=abs(A.geno-2), B.geno=abs(B.geno-2))]
    tmp.ag[,fRR:=nRR/n]
    tmp.ag[,fRA:=nRA/n]
    tmp.ag[,fAA:=nAA/n]

    tmp.ag[,sample.id:=seqGetData(genofile, "sample.id")[i]]
    tmp.ag
  }

  f1.set <- het.count[n>1000][A.geno==0 & B.geno==2 & fRA>.9]$sample.id


  het.count.ag <- het.count[,list(RR=c(mean(fRR[A.geno==0]), mean(fRR[B.geno==0])),
                                  RA=c(mean(fRA[A.geno==1]), mean(fRA[B.geno==1])),
                                  AA=c(mean(fAA[A.geno==2]), mean(fAA[B.geno==2])),
                                  n=sum(n),
                                  SC=c("A", "B")),
                             list(sample.id)]
  het.count.ag[,mu:=(RR+RA+AA)/3]


### make small, dummy vcf
    seqResetFilter(genofile)

    seqSetFilter(genofile, sample.id=c(het.count.ag[SC=="A"][which.max(mu)]$sample.id,
                                       het.count.ag[SC=="B"][which.max(mu)]$sample.id,
                                       f1.set[1:10]),
                           variant.id=snp.dt$id)

     seqGDS2VCF(genofile, vcf.fn="/scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.10f1.vcf.gz", info.var=character(0), fmt.var=character(0))

     ped <- data.table(fam="family1",
                       iid=f1.set[1:10],
                       pid=het.count.ag[SC=="B"][which.max(mu)]$sample.id,
                       mid=het.count.ag[SC=="A"][which.max(mu)]$sample.id,
                       foo="N", bar="A")

     write.table(ped, sep="\t", quote=F, row.names=F, col.names=F, file="/scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.10f1.ped")

### make ped file
    ped <- data.table(fam="family1",
                      iid=f1.set,
                      pid=sc[SC=="A"]$clone[1],
                      mid=sc[SC=="B"]$clone[1],
                      foo="N", bar="A")

    write.table(ped, sep="\t", quote=F, row.names=F, col.names=F, file="/scratch/aob2x/daphnia_hwe_sims/trioPhase/allTrio.ped")





  save(het.count, file="/nv/vol186/bergland-lab/alan/f1q.Rdata")


  library(data.table)
  library(ggplot2)

  load("/mnt/sammas_storage/bergland-lab/alan/f1q.Rdata")

  hist(het.count[A.geno==0 & B.geno==2]$fRA)
  abline(v=.90, col="red")


  het.count[sample.id%in%het.count[n>10][A.geno==0 & B.geno==2 & fRA>.9]$sample.id]
