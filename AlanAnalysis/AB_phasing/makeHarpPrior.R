#ijob -c1 -p standard -A berglandlab
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R


library(data.table)
library(foreach)

dat <- fread("/scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.csv")
setnames(dat, c(1,2), c("chr", "pos"))

dat <- dat[!(A=="1/1" & B=="1/1")][!(A=="0/0" & B=="0/0")]

dat.p <- dat[,list(Ref=REF,
                    A1=c(ALT, REF)[as.numeric(substr(A, 0, 1))+1],
                    A2=c(ALT, REF)[as.numeric(substr(A, 3, 3))+1],
                    B1=c(ALT, REF)[as.numeric(substr(B, 0, 1))+1],
                    B2=c(ALT, REF)[as.numeric(substr(B, 3, 3))+1],
                    Coverage=4),
                list(chr, pos)]

dat.ag <- dat[,list(.N, max=max(pos)), chr]
dat.ag <- dat.ag[N>1000]
dat.ag[,id:=1:12]

setkey(dat.p, chr)
foreach(i= dat.ag$chr)%do%{
  print(i)
  # i <- dat.ag[N>1000]$chr[1]
  tmp <- dat.p[J(i)][,-"chr",with=F]

  setnames(tmp, "pos", i)

  write.table(tmp, file=paste("/scratch/aob2x/daphnia_hwe_sims/harp_pools/priors/", i, ".csv", sep=""), quote=F, row.names=F, sep=",")
  system(paste("sed 's/$/,/g' /scratch/aob2x/daphnia_hwe_sims/harp_pools/priors/", i, ".csv | gzip -c - > /scratch/aob2x/daphnia_hwe_sims/harp_pools/priors/", i, ".csv.gz", sep=""))
}


### make job id file
pools <- system("ls /scratch/aob2x/daphnia_hwe_sims/harp_pools/bams/*bam", intern=T)
jobs <- foreach(p=pools, .combine="rbind")%do%{
  tmp <- dat.ag[,c("chr", "max"),with=F]
  tmp[,bam:=p]
}
jobs[,id:=c(1:dim(jobs)[1])]
jobs <- jobs[,c("id", "chr", "max", "bam"),with=F]


write.table(jobs, file="/scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId", quote=F, row.names=F, col.names=T)
