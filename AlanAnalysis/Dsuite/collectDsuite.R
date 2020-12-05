
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)


### f stats
  f.shapeit <- fread("/scratch/aob2x/daphnia_hwe_sims/dsuite/shapeit.bbaa.txt")
  f.shapeit[,data:="shapeit"]

  f.orig <- fread("/scratch/aob2x/daphnia_hwe_sims/dsuite/orig.bbaa.txt")
  f.orig[,data:="orig"]

  f <- rbind(f.shapeit, f.orig)
  setnames(f, "Z-score", "z")
  fw <- dcast(f, P1 + P2 + P3 ~ data, value.var="z")

  save(f, fw, file="~/fstat.Rdata")




### scp aob2x@rivanna.hpc.virginia.edu:~/fstat.Rdata ~/.

  library(data.table)
  library(ggplot2)

  load("~/fstat.Rdata")

  ggplot(data=fw, aes(x=orig, y=shapeit, color=as.factor(P3=="pulicaria"))) + geom_point() + geom_abline(slope=1, intercept=0)
