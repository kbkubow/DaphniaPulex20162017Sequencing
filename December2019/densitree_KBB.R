#!/usr/bin/env Rscript

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(phangorn)
  library(doMC)
  registerDoMC(20)

### open genotype file
  genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds",
                      allow.duplicate=TRUE)

### Load LD pruned SNP file to use for Filtering
  filtsnptb <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexandObtusaC/finalsetsnpset01pulex_table_20200207")
  colnames(filtsnptb) <- c("oldvariantids", "chr", "pos", "olddp")

### Filter SNPs
  snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
      chr = seqGetData(genofile, "chromosome"),
      pos = seqGetData(genofile, "position"))
  setkey(snps, chr, pos)
  setkey(filtsnptb, chr, pos)
  msnps <- merge(filtsnptb, snps)
  msnps[,final.use:=T]

### Load superclone file (does not contain Pulicaria)
  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexandObtusaC/CloneInfoFilePulexandObtusa_withmedrd_20200207")

### Do some filtering on read depth and independence
  scrd5 <- sc[medrd>14 & Nonindependent==0]
  #scrd5 <- sc[medrd>4 & Nonindependent==0 & Species=="pulex"]
  subscrd5 <- scrd5[, c("clone", "SC", "population", "year", "Sex", "Species"), with=FALSE]
  #clonestouse <- subscrd5

### Optional subsetting of superclones
  others <- subscrd5[SC=="OO"]
  scs <- subscrd5[SC!="OO"]
  listofsc <- unique(scs$SC)
  setkey(scs, SC)


  scsub <- foreach(sc.i=listofsc, .combine="rbind")%dopar%{
    #sc.i <- "A"
    sctmp <- scs[J(sc.i)]
    sctmpsub <- if(dim(sctmp[1] > 2, sctmp[sample(nrow(sctmp), 3),]
    sctmpsub
    }

    data.table(chr=chr.i, start=seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size),
                          stop =seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size) + step.size)
  }


### Add in pulicaria
  pulicaria <- data.table(clone=c("2018_Pulicaria_Pond21_22",
      "2018_Pulicaria_Pond22_21", "2018_Pulicaria_Pond22_53",
      "2018_Pulicaria_Pond22_62", "2018_Pulicaria_Pond22_72"), SC=c("ZZ"), population=c("Pond21", "Pond22",
      "Pond22", "Pond22", "Pond22"), year=c("2018"), Sex=c("female"), Species=c("pulicaria"))
  clonestouse <- rbind(subscrd5, pulicaria)
  clonestouseids <- clonestouse$clone

### figure out who is root
  setkey(clonestouse, clone)
  sc.root <- data.table(clone=(seqGetData(genofile, "sample.id")))
  sc.root

  #use.root <- "April_2017_Dbarb_11"
  use.root <- "Spring_2016_W1_1.1"

### make windows
  setkey(msnps, chr, pos)

  chr.ag <- msnps[,list(start = min(pos), stop=max(pos)), list(chr)]
  chr.ag$length <- chr.ag$stop-chr.ag$start
  chrtouse <- chr.ag$chr[chr.ag$length>1000000]


  #window.size <- 1000000
  #step.size <- 1000000
  window.size <- 2500000
  step.size <- 2500000
  wins <- foreach(chr.i=chrtouse, .combine="rbind")%dopar%{
    #chr.i <- "Scaffold_1863_HRSCAF_2081"
    snp.dt.ag <- msnps[J(chr.i)][(final.use), list(start=min(pos), stop=max(pos))]

    data.table(chr=chr.i, start=seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size),
                          stop =seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size) + step.size)
  }

  ### Maybe look at how many SNPs go into each window/tree? Perhaps there are some low ones that are throwing things off?
  msnpcount <- foreach(win.i = c(1:dim(wins)[1]), .combine="rbind", .errorhandling="remove")%dopar%{

    print(paste("Getting data: ", win.i, sep=""))

    snp.dt.tmp <- msnps[J(wins[win.i]$chr)][pos>=wins[win.i]$start & pos<=wins[win.i]$stop][(final.use)]

    tmp <- data.table(window=win.i, numsnps=dim(snp.dt.tmp)[1])
    tmp

  }

  wins$numsnps <- msnpcount$numsnps
  winstouse <- wins[numsnps>800]


### run windows
  m <- foreach(win.i = c(1:dim(winstouse)[1]), .errorhandling="remove")%dopar%{
  #m <- foreach(win.i =1:20, .errorhandling="remove")%dopar%{

    print(paste("Getting data: ", win.i, sep=""))

    snp.dt.tmp <- msnps[J(winstouse[win.i]$chr)][pos>=winstouse[win.i]$start & pos<=winstouse[win.i]$stop][(final.use)]

    seqResetFilter(genofile)
    seqSetFilter(genofile, variant.id=snp.dt.tmp$variant.ids, sample.id=clonestouse[!is.na(population)]$clone)

    dosage <- seqGetData(genofile, "$dosage")

    rownames(dosage) <- seqGetData(genofile, "sample.id")

    print(paste("Making Tree: ", win.i, sep=""))
    nj.o  <- njs(dist(dosage))

    nj.o <- root(phy=nj.o, outgroup=use.root)

  }

  class(m) <- "multiPhylo"

  setkey(clonestouse, clone)
  clones <- clonestouse[J(m[[1]]$tip.label)]
  #save(m, clones, file="/njofilt.Rdata")

### plot
  #load("njo.Rdata")
  clones.mat <- as.matrix(clones[,"population", with=T])
  rownames(clones.mat) <- clones$clone

  clones[population=="D8", col:="blue"]
  clones[population=="DBunk", col:="red"]
  clones[population!="D8" & population!="DBunk", col:="black"]

  clones[,u:=as.numeric(as.factor(clone))]

  ### clonal labeling
    clones[SC=="A", lab:="  A"]
    clones[SC=="C", lab:="    C"]
    clones[SC=="B", lab:="      B"]
    clones[SC=="D", lab:="        D"]
    clones[SC=="G", lab:="          G"]
    clones[is.na(lab), lab:=""]
    setkey(clones,  clone)

  ### pond labeling
    clones[,pondf := as.numeric(as.factor(population))]
    clones[,lab := apply(cbind(clones$population, clones$pondf, clones$u, lab), 1, function(x) paste(x[4],
                                              paste(rep("       ", as.numeric(x[2])), collapse=""),
                                                                        x[1], "  ", x[3], sep=""))]

  m.new <- foreach(m.i=1:length(m))%do%{
    tmp <-  m[[m.i]]
    tmp$tip.label <- clones[J(tmp$tip.label)]$lab

    tmp <- root(tmp, clones[clone==use.root]$lab)
    tmp
  }
  class(m.new) <- "multiPhylo"

  plot(m.new[[1]])
  plot(root(m.new[[1]], 1))

  pdf("~/densiTree.pdf", height=5, w=9)
  densiTree(m.new, type="cladogram", scaleX=T, use.edge=F, root.edge=T, x.lim=c(400, 0))
  dev.off()
