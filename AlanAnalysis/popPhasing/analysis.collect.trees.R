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
  cdl.o <- foreach(i=fn[1:5])%dopar%{
    message(i)
    #i <- fn[1]
    load(i)
    #return(cdl)
    cdl[,window:=tstrsplit(i1, ";")[[2]]]

    cdl[,cd_bin:=round(cd, 5)]
    #cdl[,cd_bin:=factor(cd_bin, seq(from=0, to=.1, by=.001))]
    cdl[,list(n=.N), list(sp.group, pond.group, window, cd_bin)]

  }

  cdl.o <- rbindlist(cdl.o)
  cdl.genome <- cdl.o[,list(n=sum(n)), list(sp.group, pond.group, cd_bin)][!is.na(sp.group) & !is.na(pond.group)]

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
