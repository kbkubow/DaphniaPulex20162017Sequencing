#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### makes rabbit input

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone=
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### which F1s?
  f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.delim")

### open GDS
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")

### make snp.dt
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       key="chr")

### make large input file
  chr.i="Scaffold_9199_HRSCAF_10755"
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

  genomat <- genomat[sample(c(1:dim(genomat)[1]), 1000)]
  genomat <- genomat[order(id)]

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
  pos <- matrix(c("chromosome", seq(from=0, to=100, length.out=dim(genomat)[1])), nrow=1)

  header <- do.call("rbind", list(marker, chr, pos))

  out <- do.call("rbind", list(header, parents, offspring))



  out.fn <- paste("/scratch/aob2x/", chr.i, ".all.in", sep="")

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






  genomat[ind=="0", ind:="2N"] ### I think this is wrong shoudl be 1N, and other 2N
  genomat[ind=="1", ind:=sample(c("1N","2N"), length(ind), replace=T)]
  genomat[ind=="2", ind:= "1N"]
  genomat[is.na(ind), ind:= "NN"]












  getOutput <- function(ind.i, chr.i) {

    #ind.i="March20_2018_D8_18030"; chr.i="Scaffold_9199_HRSCAF_10755"
    seqSetFilter(genofile,
                sample.id=c(sc[SC=="A"][which.max(medrd)]$clone,
                            sc[SC=="C"][which.max(medrd)]$clone,
                            ind.i),
                variant.id=snp.dt[J(chr.i)][numAlleles==2]$id)

    genomat <- as.data.table(t(seqGetData(genofile, "$dosage")))
    setnames(genomat, seqGetData(genofile, "sample.id"))

    setnames(genomat, sc[SC=="A"][which.max(medrd)]$clone, "A")
    setnames(genomat, sc[SC=="C"][which.max(medrd)]$clone, "C")
    setnames(genomat, ind.i, "ind")

    genomat[,id:=as.character(seqGetData(genofile, "variant.id"))]
    genomat <- merge(genomat, vcf[,c("id", "A", "C", "Acode", "Ccode", "cm", "chr"), with=F], by="id")

    #table(genomat$C.x, genomat$C.y)
    genomat[,ind.orig:=ind]
    genomat[,ind:=as.character(ind)]

    genomat[ind=="0", ind:="2N"] ### I think this is wrong shoudl be 1N, and other 2N
    genomat[ind=="1", ind:=sample(c("1N","2N"), length(ind), replace=T)]
    genomat[ind=="2", ind:= "1N"]
    genomat[is.na(ind), ind:= "NN"]

    table(genomat$ind, genomat$ind.orig)
    table(genomat$Acode, genomat$A.x)

    table(genomat$Acode, genomat$Ccode, genomat$ind)
    table(genomat$C, genomat$A, genomat$ind)


  ### make header info

    marker <- matrix(c("marker", genomat$id), nrow=1)
    chr <- matrix(c("chromosome", as.numeric(as.factor(genomat$chr))), nrow=1)
    pos <- matrix(c("pos(cM)", genomat$cm), nrow=1)

    header <- do.call("rbind", list(marker, chr, pos))

    data <- cbind(matrix(c("C", "A", ind.i), ncol=1), t(as.matrix(genomat[,c("Ccode", "Acode", "ind")])))

    gmh <- rbind(header, data)


  ### write output, per chromosome per individual
    #out.fn <- paste("/scratch/aob2x/daphnia_hwe_sims/trioPhase/rabbitIn/", chr.i, ".", ind.i, ".in", sep="")

    out.fn <- paste("/scratch/aob2x/", chr.i, ".", ind.i, ".in", sep="")

    writeLines( paste("#founders,",2, sep=""),
                 con=out.fn
               )
    options(scipen=999)

     write.table(gmh,
                 file=out.fn,
                 quote=FALSE,
                 row.names=FALSE,
                 col.names=FALSE,
                 sep=",",
                 na="N",
                 append=TRUE)
   }


   foreach(chr.i=unique(snp.dt$chr))%do%{
     foreach(ind.i=f1s$cloneid)%do%{
       message(paste(chr.i, ind.i, sep=" / "))
       getOutput(ind.i, chr.i)
     }
   }
