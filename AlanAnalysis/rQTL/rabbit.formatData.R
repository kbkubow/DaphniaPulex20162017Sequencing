#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### makes rabbit input
  args = commandArgs(trailingOnly=TRUE)
  chr.i <- as.character(args[1])
  maxcM <- as.numeric(args[2])
  #chr.i <- "Scaffold_9199_HRSCAF_10755"

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone=
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### which F1s?
  #f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.onlyPheno.delim")
  f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.allF1s.delim")

### open GDS
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

### make snp.dt
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       key="chr")

  #snp.dt.ag <- snp.dt[,.N,chr]
  #write.table(snp.dt.ag, file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase/chrs.csv", quote=F, row.names=T, col.names=F, sep=",")

### make large input file

  seqSetFilter(genofile,
              sample.id=c(sc[SC=="A"][which.max(medrd)]$clone,
                          sc[SC=="C"][which.max(medrd)]$clone,
                          f1s$cloneid),
              variant.id=snp.dt[J(chr.i)][numAlleles==2]$id)

  genomat <- as.data.table(t(seqGetData(genofile, "$dosage")))
  setnames(genomat, seqGetData(genofile, "sample.id"))

  setnames(genomat, sc[SC=="A"][which.max(medrd)]$clone, "A")
  setnames(genomat, sc[SC=="C"][which.max(medrd)]$clone, "C")

  genomat[,id:=seqGetData(genofile, "variant.id")]

  genomat <- genomat[(A!=0 & C!=2) | (A!=2 & C!=0)]

  ### random sample
    #genomat <- genomat[sample(c(1:dim(genomat)[1]), 5000)]
    #genomat <- genomat[order(id)]

  ### most informative
    genomat <- genomat[(A==1 & C==0) | (A==1 & C==2) | (A==0 & C==1) | (A==2 & C==1)]
    setkey(genomat, id)
    genomat <- genomat[J(sample(genomat$id, 5000))]
    genomat <- genomat[order(id)]

    #table(genomat$A, genomat$C)

  parents <- foreach(ind.i=c("A", "C"), .combine="rbind")%do%{
    tmp <- t(as.matrix(genomat[,ind.i, with=F]))
    tmp[tmp=="2"] <- "22"
    tmp[tmp=="1"] <- "12"
    tmp[tmp=="0"] <- "11"

    cbind(matrix(ind.i, ncol=1), tmp)
  }

  offspring <- foreach(ind.i=f1s$cloneid, .combine="rbind", .errorhandling="remove")%do%{
    tmp <- t(as.matrix(genomat[,ind.i, with=F]))
    tmp[tmp=="2"] <- "2N"
    tmp[tmp=="1"] <- sample(c("1N","2N"), dim(tmp)[1], replace=T)
    tmp[tmp=="0"] <- "1N"
    tmp[is.na(tmp)] <- "NN"
    cbind(matrix(ind.i, ncol=1), tmp)
  }

  marker <- matrix(c("marker", genomat$id), nrow=1)
  #chr <- matrix(c("chromosome", rep(NA, dim(genomat)[1])), nrow=1)
  #pos <- matrix(c("pos(cM)", rep(NA, dim(genomat)[1])), nrow=1)
  chr <- matrix(c("chromosome", rep(as.numeric(as.factor(chr.i)), dim(genomat)[1])), nrow=1)
  pos <- matrix(c("chromosome", seq(from=0, to=maxcM, length.out=dim(genomat)[1])), nrow=1)

  header <- do.call("rbind", list(marker, chr, pos))

  out <- do.call("rbind", list(header, parents, offspring))



  out.fn <- paste("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_", maxcM, "cm/", chr.i, ".all.in", sep="")

  writeLines( paste("#founders,",2, sep=""),
               con=out.fn
             )
  options(scipen=999)

   write.table(out,
               file=out.fn,
               quote=FALSE,
               row.names=FALSE,
               col.names=FALSE,
               sep=",",
               na="NA",
               append=TRUE)
