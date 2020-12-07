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

  cdl.list <- foreach(i=fn)%dopar%{
    message(i)
    #i <- fn[1]
    load(i)
    #return(cdl)
    cdl[,window:=tstrsplit(i1, ";")[[2]]]

    cdl[,cd_bin:=round(cd, 5)]
    #cdl[,cd_bin:=factor(cd_bin, seq(from=0, to=.1, by=.001))]

    njo$tip.label <- gsub(".fa", "", tstrsplit(njo$tip.label, ";")[[1]])

    njo <- root(njo, njo$tip.label[grepl("obtusa", njo$tip.label)])

    o.tmp <- list(cdl[,list(n=.N), list(sp.group, pond.group, window, cd_bin)],
                  njo)

    names(o.tmp) <- c(cdl[1]$window, cdl[1]$window)

    o.tmp
  }

  cdl.o <- lapply(cdl.list, function(x) x[[1]])
  cdl.o <- rbindlist(cdl.o)
  cdl.genome <- cdl.o[,list(n=sum(n)), list(sp.group, pond.group, cd_bin)][!is.na(sp.group) & !is.na(pond.group)]


  cdl.tree <- lapply(cdl.list, function(x) x[[2]])
  names(cdl.tree) <- lapply(cdl.list, function(x) names(x)[2])


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


  save(cdl.genome, cdl.qtl, cdl.o, cdl.tree, cdl.list, file="~/cdlo.Rdata")
