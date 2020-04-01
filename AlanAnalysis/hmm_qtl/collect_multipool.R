### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(doMC)

### import
  fl <- system("ls /scratch/aob2x/daphnia_hwe_sims/multipool/output/*.out", intern=T)

  o <- foreach(i=fl)%do%{
    #i <- fl[1]
    tmp <- fread(i)
    tmp[,chr:=gsub(".out", "", last(tstrsplit(i, "/")))]
  }
  o <- rbindlist(o)
  setnames(o, names(o), c("bin", "mle", "lod", "chr"))

### export
  save(o, file="/nv/vol186/bergland-lab/alan/multipool.Rdata")

### import

  library(data.table)
  library(ggplot2)
  library(foreach)
  library(doMC)

  load("/mnt/sammas_storage/bergland-lab/alan/multipool.Rdata")

  ggplot(data=o, aes(x=bin, y=lod, color=chr)) + geom_line() + facet_wrap(~chr, nrow=1)
