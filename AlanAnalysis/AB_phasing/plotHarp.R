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
  save(o, file="/scratch/aob2x/daphnia_hwe_sims/harp_pools/summarizedOut/harpWins_250K.Rdata")
  save(o, file="/nv/vol186/bergland-lab/alan/harpWins_250K.Rdata")


### libraries
  library(ggplot2)
  library(data.table)

  load("/mnt/sammas_storage/bergland-lab/alan/harpWins.Rdata")
  o[,mid:=(start+stop)/2]

  o[,mid:=(start+stop)/2]
  o.ag <- o[,.N,list(chr,mid)]
  setkey(o.ag, chr, mid)
  o.ag[,i:=1:dim(o.ag)[1]]

  setkey(o, chr, mid)
  o <- merge(o, o.ag[,-"N",with=F])

  setkey(o, chr, mid)

  ol <- melt(o, id.vars=c("chr", "start", "stop", "pool", "mid", "i"), variable.name = "haplotype", value.name = "f")

  ol.ag <- ol[,list(f.delta=(qlogis(f)-qlogis(mean(f))), pool=pool), list(chr, mid, haplotype, i)]
  ol.ag.ag <- ol.ag[,list(f.delta.mu=f.delta[which.min(abs(f.delta))]), list(chr, mid, haplotype, i, pool=ifelse(grepl("Male", pool), "male", "pe"))]


  ggplot(data=ol.ag.ag, aes(x=i, y=f.delta.mu, group=pool, color=haplotype)) + facet_wrap(~haplotype, ncol=1) +
  geom_vline(xintercept=ol.ag[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(mid-8660157))]$i) +  geom_line()


  ii <- ol.ag[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(mid-8660157))]$i
  ii <- ol.ag.ag[f.delta.mu< -1]$i
  ggplot(data=ol.ag.ag[i>=(ii-10) & i<=(ii+10)], aes(x=i, y=f.delta.mu, group=pool, color=pool)) + geom_line() + facet_wrap(~haplotype, ncol=1)

wins[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(start-8660157))]$i









  o[,mid:=(start+stop)/2]
  o.ag <- o[,.N,list(chr,mid)]
  setkey(o.ag, chr, mid)
  o.ag[,i:=1:dim(o.ag)[1]]

  setkey(o, chr, mid)
  o <- merge(o, o.ag[,-"N",with=F])

  ol <- melt(o, id.vars=c("chr", "start", "stop", "pool", "mid", "i"), variable.name = "haplotype", value.name = "f")

  ggplot(data=ol, aes(y=f, x=haplotype, color=interaction(haplotype, pool))) + geom_hline(yintercept=0.25) + geom_boxplot()

  ggplot(data=ol, aes(y=f, x=i, color=chr, group=haplotype)) + geom_line() + facet_wrap(~pool, ncol=1)

  ggplot(data=ol[chr=="Scaffold_7757_HRSCAF_8726"], aes(y=f, x=i, color=haplotype, group=haplotype)) +
  geom_line() + facet_wrap(~pool, ncol=1)

  ol.ag <- ol[,list(male=mean(f[grepl("Male", pool)]), pe=mean(f[grepl("PE", pool)])), list(chr, mid, i, haplotype)]

  ol.ag[,diff:=male-pe]
  ol.ag[,nDiff:=qlogis(male) - qlogis(pe)]

  ggplot(data=ol.ag, aes(x=i, y=diff, group=haplotype, color=haplotype)) + geom_line()
  ggplot(data=ol.ag[chr=="Scaffold_7757_HRSCAF_8726"], aes(x=i, y=nDiff, group=haplotype, color=haplotype)) + geom_line()


  o.ag2 <- ol.ag[,list(male.or=(male[haplotype=="pA1"] /male[haplotype=="pA2"]) / (male[haplotype=="pB1"] /male[haplotype=="pB2"]),
                       pe.or=(pe[haplotype=="pA1"] /pe[haplotype=="pA2"]) / (pe[haplotype=="pB1"] /pe[haplotype=="pB2"])), list(chr, mid, i)]


   o.ag2 <- ol.ag[,list(male.or=(male[haplotype=="pA1"] /male[haplotype=="pB1"]) / (male[haplotype=="pA2"] /male[haplotype=="pB2"]),
                        pe.or=(pe[haplotype=="pA1"] /pe[haplotype=="pB1"]) / (pe[haplotype=="pA2"] /pe[haplotype=="pB2"])), list(chr, mid, i)]




   o.ag2 <- ol.ag[,list(a.or=(male[haplotype=="pA1"] /pe[haplotype=="pA1"]) / (male[haplotype=="pA2"] /pe[haplotype=="pA2"]),
                        b.or=(male[haplotype=="pB1"] /pe[haplotype=="pB1"]) / (male[haplotype=="pB2"] /pe[haplotype=="pB2"])), list(chr, mid, i)]


  o.ag2 <- ol.ag[,list(a.or=(male[haplotype=="pA1"] /male[haplotype=="pA2"]) / (pe[haplotype=="pA1"] /pe[haplotype=="pA2"]),
                       b.or=(male[haplotype=="pB1"] /male[haplotype=="pB2"]) / (pe[haplotype=="pB1"] /pe[haplotype=="pB2"])), list(chr, mid, i)]



  o.ag2[,a.or:=log2(a.or)]
  o.ag2[,b.or:=log2(b.or)]

     o.ag2[,diff:=(male.or)


   o.ag2[,diff:=log2(a.or/b.or)]


 ggplot(data=o.ag2[chr=="Scaffold_7757_HRSCAF_8726"], aes(x=i, y=b.or)) + geom_line()
 ggplot(data=o.ag2, aes(x=i, y=b.or, color=chr)) + geom_line()
 ggplot(data=o.ag2, aes(x=i, y=diff, color=chr)) + geom_line()



   o.ag2 <- ol.ag[,list(or=(male[haplotype=="pA1"]  + male[haplotype=="pA2"]) / (pe[haplotype=="pA1"] + pe[haplotype=="pA2"]) /
                            (male[haplotype=="pB1"] + male[haplotype=="pB2"]) / (pe[haplotype=="pB1"] + pe[haplotype=="pB2"])), list(chr, mid, i)]

  ggplot(data=o.ag2, aes(x=i, y=log2(or), color=chr)) + geom_line()



ol.ag[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(mid-8660157))]
ol.ag[i==7187]



o.ag2 <- ol.ag[,list(male.ss=sum((male-.25)^2), pe.ss=sum((pe-.25)^2)), list(chr, mid, i)]
ggplot(data=o.ag2, aes(x=i, y=male.ss, color=chr)) + geom_line()

ggplot(ol.ag, aes(x=male, y=pe, color=haplotype)) + geom_hex() + facet_wrap(~haplotype)


ol.ag.ag <- ol.ag[,list(diff=)]
