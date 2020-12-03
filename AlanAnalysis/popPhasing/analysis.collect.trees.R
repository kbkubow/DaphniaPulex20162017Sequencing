#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


### libraries
  library(data.table)
  library(foreach)
  library(ape)
  library(doMC)
  registerDoMC(20)

### input files: genome-wide distribution


  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/popPhase/trees/", "Rdata", full.names=T)
  length(fn)
  cdl.o <- foreach(i=fn)%dopar%{
    message(i)
    #i <- fn[1]
    load(i)
    #return(cdl)
    cdl[,window:=tstrsplit(i1, ";")[[2]]]

    cdl[,list(cd_mean=mean(cd), cd_sd=sd(cd)), list(sp.group, pond.group, window)]
  }

  cdl.o <- rbindlist(cdl.o)
  cdl.genome <- cdl.o[,list(cd_mean=mean(cd), cd_sd=sd(cd)), list(sp.group, pond.group, window)]

### qtl
  wins <- fread("/scratch/aob2x/daphnia_hwe_sims/popPhase/windows.delim")

  target.fn <- sapply(tail(wins, 14)$V2, function(x) fn[grepl(x, fn)])

  length(target.fn)
  cdl.qtl <- foreach(i=target.fn)%dopar%{
    message(i)
    #i <- fn[1]
    load(i)
    #return(cdl)
    cdl[,window:=tstrsplit(i1, ";")[[2]]]

    #cdl[,list(cd_mean=mean(cd), cd_sd=sd(cd)), list(sp.group, pond.group, window)]
  }
  cdl.qtl <- rbindlist(cdl.qtl)


### save


  save(cdl.genome, cdl.qtl, cdl.o, file="~/cdlo.Rdata")
