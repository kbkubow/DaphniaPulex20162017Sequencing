#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(SeqArray)

### load parents and format parental haplotypes
  vcf <- fread("/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.noNA.vcf", skip="#CHROM", na.strings="./.")
  vcf[,A1:=as.numeric(tstrsplit(A, "\\||\\/")[[1]]) + 1]
  vcf[,A2:=as.numeric(tstrsplit(A, "\\||\\/")[[2]]) + 1]
  vcf[,C1:=as.numeric(tstrsplit(C, "\\||\\/")[[1]]) + 1]
  vcf[,C2:=as.numeric(tstrsplit(C, "\\||\\/")[[2]]) + 1]
  vcf[,id:=tstrsplit(ID, "_")[[2]]]
  setnames(vcf, "#CHROM", "chr")

  vcf.ag <- vcf[,list(maxPos=max(POS)), list(chr)]
  setkey(vcf, chr)
  setkey(vcf.ag, chr)

  vcf <- merge(vcf, vcf.ag)
  vcf[,cm:=(POS/maxPos) * 100]

### which F1s?
  f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.delim")

### open GDS
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")

  getOutput <- function(ind.i, chr.i) {

    #ind.i="AxB_R4_P17_B"; chr.i="Scaffold_2373_HRSCAF_2879"
    seqSetFilter(genofile, sample.id=ind.i, variant.id=vcf[chr==chr.i]$id)

    genomat <- seqGetData(genofile, "$dosage_alt")
    genomat[is.na(genomat)] <- "NN"
    genomat[genomat=="0"] <- "11"
    genomat[genomat=="1"] <- "12"
    genomat[genomat=="2"] <- "22"

    f1.genomat <- cbind(matrix(seqGetData(genofile, "sample.id"), ncol=1), genomat)

  ### merge parents and offspring
    parent.genomat <- cbind(matrix(c("A1", "A2", "C1", "C2"), ncol=1),
                            t(as.matrix(vcf[chr==chr.i,c("A1", "A2", "C1", "C2"), with=F])))

    gm <- rbind(parent.genomat, f1.genomat)

  ### make header info



    marker <- matrix(c("marker", vcf[chr==chr.i]$id), nrow=1)
    chr <- matrix(c("chromosome", as.numeric(as.factor(vcf[chr==chr.i]$chr))), nrow=1)
    pos <- matrix(c("pos(cM)", vcf[chr==chr.i]$cm), nrow=1)

    header <- do.call("rbind", list(marker, chr, pos))

    gmh <- rbind(header, gm)


  ### write output, per chromosome per individual
    out.fn <- paste("/scratch/aob2x/daphnia_hwe_sims/trioPhase/rabbitIn/", chr.i, ".", ind.i, ".in", sep="")

    writeLines( paste("#founders,",4, sep=""),
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
