#ijob -c1 -p standard -A berglandlab
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)
  library(foreach)

### load harp output
  dirs <- system("ls /scratch/aob2x/daphnia_hwe_sims/harp_pools/out/", intern=T)

  o <- foreach(i=dirs, .combine="rbind")%do%{
    #i <- dirs[1]
    print(i)
    tmp <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/harp_pools/out/", i, "/", i, ".freqs", sep=""))
    setnames(tmp, names(tmp), c("chr", "start", "stop", "pA1", "pA2", "pB1", "pB2"))
    tmp[,pool:=tstrsplit(i, "\\.")[[1]]]

    tmp
  }

### save
  save(o, file="/scratch/aob2x/daphnia_hwe_sims/harp_pools/summarizedOut/harpWins.Rdata")
  save(o, file="/nv/vol186/bergland-lab/alan/harpWins.Rdata")


### libraries
  library(ggplot2)
  library(data.table)

  load("/mnt/sammas_storage/bergland-lab/alan/harpWins.Rdata")
  o[,mid:=(start+stop)/2]
  o.ag <- o[,.N,list(chr,mid)]
  setkey(o.ag, chr, mid)
  o.ag[,i:=1:dim(o.ag)[1]]

  setkey(o, chr, mid)
  o <- merge(o, o.ag[,-"N",with=F])

  ol <- melt(o, id.vars=c("chr", "start", "stop", "pool", "mid", "i"), variable.name = "haplotype", value.name = "f")

  ggplot(data=ol, aes(y=f, x=haplotype, color=haplotype)) + geom_hline(yintercept=0.25) + geom_boxplot()

  ggplot(data=ol, aes(y=f, x=i, color=chr, group=haplotype)) + geom_line() + facet_wrap(~pool, ncol=1)

  ggplot(data=ol[chr=="Scaffold_7757_HRSCAF_8726"], aes(y=f, x=i, color=haplotype, group=haplotype)) +
  geom_line() + facet_wrap(~pool, ncol=1)
