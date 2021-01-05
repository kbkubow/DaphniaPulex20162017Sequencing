
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)


### f stats
  f.shapeit <- fread("/scratch/aob2x/daphnia_hwe_sims/dsuite/shapeit.bbaa.txt")
  f.shapeit[,data:="shapeit"]

  f.orig <- fread("/scratch/aob2x/daphnia_hwe_sims/dsuite/orig.bbaa.txt")
  f.orig[,data:="orig"]

  f.hybrid <- fread("/scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy_BBAA.txt")
  f.hybrid[,data:="hybrid"]
    #setnames(f.hybrid, "Z-score", "z")
    #f.hybrid[P3=="pulicaria"][order(z)][,c('P1', 'P2', 'P3', 'z', 'Dstatistic', 'BBAA', 'ABBA', 'BABA'), with=F]

  f.hybrid.sp3 <- fread("/scratch/aob2x/daphnia_hwe_sims/dsuite/hybrid_strategy.3species_BBAA.txt")
  f.hybrid.sp3[,data:="hybrid.sp3"]

  f <- rbindlist(list(f.shapeit, f.orig, f.hybrid, f.hybrid.sp3))

  setnames(f, "Z-score", "z")
  fw <- dcast(f, P1 + P2 + P3 ~ data, value.var="z")

  save(f, fw, file="~/fstat.Rdata")




### scp aob2x@rivanna.hpc.virginia.edu:~/fstat.Rdata ~/.

  library(data.table)
  library(ggplot2)

  load("~/fstat.Rdata")

  ggplot(data=fw, aes(x=orig, y=shapeit, color=as.factor(P3=="pulicaria"))) + geom_point() + geom_abline(slope=1, intercept=0)

  ggplot(data=fw, aes(x=shapeit, y=hybrid.sp3, color=as.factor(P3=="pulicaria"))) + geom_point() + geom_abline(slope=1, intercept=0)

  setnames(f, "p-value", "p")
  setnames(f, "f4-ratio", "f4")

  f[data=="hybrid.sp3", pa:=p.adjust(p)]
  f[data=="hybrid.sp3"][P3=="pulicaria"][order(z)][,c('P1', 'P2', 'P3', 'Dstatistic', 'pa', 'f4-ratio', 'ABBA', 'BABA'), with=F]
 f[data=="hybrid.sp3"][f4>0.2]


 ggplot(data=f[data=="hybrid.sp3"], aes(x=-log10(p), y=f4, color=as.factor(P3=="pulicaria"))) + geom_point()
 ggplot(data=f[data=="hybrid.sp3"], aes(x=Dstatistic, y=qlogis(f4), color=as.factor(P3=="pulicaria"))) + geom_point()
