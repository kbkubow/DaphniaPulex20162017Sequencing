### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(phangorn)
  library(doMC)
  registerDoMC(20)

### open genotype file
  #genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds",
  #                    allow.duplicate=TRUE)
  genofile <- seqOpen("/mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.biallelic.phased.gds",
                      allow.duplicate=TRUE)

### load precomputed file (HWE_simulations.prep.rivanna.R makes this file)
  load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
  sc[,year:=tstrsplit(clone, "_")[[2]]]
  sc[,pond:=gsub("Dcat", "DCat", pond)]

### if using the phased data only:
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"))
  setkey(snp.dt, id)
  snp.dt[,final.use:=T]

### figure out who is root
  setkey(sc, clone)
  sc.root <- sc[J(seqGetData(genofile, "sample.id"))]
  use.root <- sc.root[pond%in%c("W1")]$clone

### make windows
  setkey(snp.dt, chr)

  window.size <- 2000000
  step.size <- 2000000
  wins <- foreach(chr.i=unique(snp.dt[(final.use)]$chr), .combine="rbind")%dopar%{
    #chr.i <- "Scaffold_1863_HRSCAF_2081"
    snp.dt.ag <- snp.dt[J(chr.i)][(final.use), list(start=min(pos), stop=max(pos))]

    data.table(chr=chr.i, start=seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size),
                          stop =seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size) + step.size)
  }

### run windows
  m <- foreach(win.i = c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{

    print(paste("Getting data: ", win.i, sep=""))

    snp.dt.tmp <- snp.dt[J(wins[win.i]$chr)][pos>=wins[win.i]$start & pos<=wins[win.i]$stop][(final.use)]

    seqResetFilter(genofile)
    seqSetFilter(genofile, variant.id=snp.dt.tmp$id, sample.id=sc.root[!is.na(pond)]$clone)

    dosage <- seqGetData(genofile, "$dosage")

    rownames(dosage) <- seqGetData(genofile, "sample.id")

    print(paste("Making Tree: ", win.i, sep=""))
    nj.o  <- nj(dist(dosage))

    nj.o <- root(phy=nj.o, outgroup=use.root)

  }
  class(m) <- "multiPhylo"

  clones <- sc[J(m[[1]]$tip.label)]

### consensus tree
  seqResetFilter(genofile)
  seqSetFilter(genofile, sample.id=sc.root[!is.na(pond)]$clone)
  cons.tree <- nj(dist(seqGetData(genofile, "$dosage")))

  save(m, clones, file="~/njo.Rdata")

### plot
  load("njo.Rdata")
  clones.mat <- as.matrix(clones[,"pond", with=T])
  rownames(clones.mat) <- clones$clone

  clones[pond=="D8", col:="blue"]
  clones[pond=="DBunk", col:="red"]
  clones[pond!="D8" & pond!="DBunk", col:="black"]

  clones[,u:=as.numeric(as.factor(clone))]

  ### clonal labeling
    clones[SC=="A", lab:="  A"]
    clones[SC=="B", lab:="    B"]
    clones[is.na(lab), lab:=""]
    setkey(clones,  clone)

  ### pond labeling
    clones[,pondf := as.numeric(as.factor(pond))]
    clones[,lab := apply(cbind(clones$pond, clones$pondf, clones$u, lab), 1, function(x) paste(x[4],
                                                                                                paste(rep("       ", as.numeric(x[2])), collapse=""),
                                                                                                x[1], "  ", x[3], sep=""))]
  ### tip labeling
    m.new <- foreach(m.i=1:length(m))%do%{
      tmp <-  m[[m.i]]
      tmp$tip.label <- clones[J(tmp$tip.label)]$lab

      tmp
    }
    class(m.new) <- "multiPhylo"


### plot
  pdf("~/densiTree.pdf", h=5, w=9)
  densiTree(m.new, type="cladogram", consensus=m.new[[1]], scaleX=T, use.edge=F, root.edge=T, x.lim=c(400, 0))
  dev.off()
