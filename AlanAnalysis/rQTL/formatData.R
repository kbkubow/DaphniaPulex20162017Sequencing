#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### makes rabbit input

### libraries
  library(data.table)
  library(SeqArray)

### load parents and format parental haplotypes
  vcf <- fread("/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.phase.noNA.vcf", skip="#CHROM", na.strings="./.")
  vcf[,A:=tstrsplit(A, ":")[[1]]]
  vcf[,C:=tstrsplit(C, ":")[[1]]]

  vcf[,A1:=as.numeric(tstrsplit(A, "\\||\\/")[[1]]) + 1]
  vcf[,A2:=as.numeric(tstrsplit(A, "\\||\\/")[[2]]) + 1]
  vcf[,C1:=as.numeric(tstrsplit(C, "\\||\\/")[[1]]) + 1]
  vcf[,C2:=as.numeric(tstrsplit(C, "\\||\\/")[[2]]) + 1]

  vcf[,Acode:=paste(A1, A2, sep="")]
  vcf[,Ccode:=paste(C1, C2, sep="")]

  vcf[,id:=tstrsplit(ID, "_")[[2]]]
  setnames(vcf, "#CHROM", "chr")

  vcf.ag <- vcf[,list(maxPos=max(POS)), list(chr)]
  setkey(vcf, chr)
  setkey(vcf.ag, chr)

  vcf <- merge(vcf, vcf.ag)
  vcf[,cm:=(POS/(maxPos+1)) * 100]



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

  getOutput <- function(ind.i, chr.i) {

    #ind.i="AxB_R4_P17_B"; chr.i="Scaffold_2373_HRSCAF_2879"
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

    genomat[ind=="0", ind:="1N"]
    genomat[ind=="1", ind:=sample(c("1N","2N"), length(ind), replace=T)]
    genomat[ind=="2", ind:= "2N"]
    genomat[is.na(ind), ind:= "NN"]

    table(genomat$ind, genomat$ind.orig)
    table(genomat$Acode, genomat$A.x)

    table(genomat$Acode, genomat$Ccode, genomat$ind)


  ### make header info

    marker <- matrix(c("marker", genomat$id), nrow=1)
    chr <- matrix(c("chromosome", as.numeric(as.factor(genomat$chr))), nrow=1)
    pos <- matrix(c("pos(cM)", genomat$cm), nrow=1)

    header <- do.call("rbind", list(marker, chr, pos))

    data <- cbind(matrix(c("A", "C", ind.i), ncol=1), t(as.matrix(genomat[,c("Acode", "Ccode", "ind")])))

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


   foreach(chr.i=unique(vcf$chr))%do%{
     foreach(ind.i=f1s$cloneid)%do%{
       message(paste(chr.i, ind.i, sep=" / "))
       getOutput(ind.i, chr.i)
     }
   }
