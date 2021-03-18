#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R
#ijob -p largemem --mem 50G -A berglandlab


### libraries
  library(data.table)
  library(foreach)
  library(ape)
  library(doMC)
  registerDoMC(15)

### input files: genome-wide distribution


  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/popPhase/trees250K/", "iqtree.Rdata", full.names=T)
  length(fn)

  cdl.list <- foreach(i=fn)%dopar%{
    #i<-fn[1]
    message(i)
    load(i)
    cdl[,chr:=paste(unlist(tstrsplit(last(tstrsplit(i, "/")), "_")[1:4]), collapse="_")]
    cdl[,start:=as.numeric(tstrsplit(window, "-")[[1]])]
    cdl[,stop:=as.numeric(tstrsplit(window, "-")[[2]])]

    list(cdl, njo)
  }

  cdl.o <- lapply(cdl.list, function(x) x[[1]])
  cdl.o <- rbindlist(cdl.o)
  cdl.genome <- cdl.o[,list(n=sum(n)), list(sp.group, pond.group, cd_bin, boot)][!is.na(sp.group) & !is.na(pond.group)]

  cdl.genome.ag <- cdl.genome[,list(cd.mean=weighted.mean(cd_bin, n)),
                               list(sp.group, pond.group, boot)][!is.na(sp.group) & !is.na(pond.group)]


  cdl.tree <- lapply(cdl.list, function(x) x[[2]])
  names(cdl.tree) <- lapply(cdl.list, function(x) names(x)[2])


### qtl
  wins <- fread("/scratch/aob2x/daphnia_hwe_sims/popPhase/windows.delim")

  target.fn <- sapply(tail(wins, 14)$V2, function(x) fn[grepl(x, fn)])

  length(target.fn)
  cdl.qtl <- foreach(i=target.fn)%dopar%{
    message(i)
    load(i)
    cdl[,chr:=paste(unlist(tstrsplit(last(tstrsplit(i, "/")), "_")[1:4]), collapse="_")]
    cdl[,start:=as.numeric(tstrsplit(window, "-")[[1]])]
    cdl[,stop:=as.numeric(tstrsplit(window, "-")[[2]])]


    cdl
  }
  cdl.qtl <- rbindlist(cdl.qtl)
  save(cdl.genome, cdl.qtl, cdl.o, cdl.tree, file="~/cdlo_250K.boot.Rdata")


### manhattan type plot
  load(file="~/cdlo_250K.boot.Rdata")
  cdl.o[,mid:=start/2 + stop/2]


  cdl.o.manhattanPlot <- cdl.o[,list(mu=sum(cd_bin*n, na.rm=T)/sum(n), min=min(cd_bin), max=max(cd_bin)),
                      list(chr, mid, sp.group, pond.group, boot)]

  cdl.o.manhattanPlot.ag <- cdl.o.manhattanPlot[,list(cd_mean=mean(mu), lci_mean=quantile(mu, .025), uci_mean=quantile(mu, .975),
                                                      cd_max=mean(max), lci_max=quantile(max, .025), uci_max=quantile(max, .975)),
                                                  list(chr, mid, sp.group, pond.group)]

  cdl.o.manhattanPlot.ag

### njo trees


  save(cdl.tree, file="~/cdlo_250K.boot.trees.Rdata")
  save(cdl.tree, file="/project/berglandlab/alan/cdlo_250K.boot.trees.Rdata")
