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

### which F1s?
  f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.delim")

### open GDS
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")

  seqSetFilter(genofile, sample.id=f1s$cloneid, variant.id=vcf$id)

  genomat <- seqGetData(genofile, "$dosage_alt")
  genomat[is.na(genomat)] <- "NN"
  genomat[genomat=="0"] <- "11"
  genomat[genomat=="1"] <- "12"
  genomat[genomat=="2"] <- "22"

  f1.genomat <- cbind(matrix(seqGetData(genofile, "sample.id"), ncol=1), genomat)

### merge parents and offspring
  parent.genomat <- cbind(matrix(c("A1", "A2", "C1", "C2"), ncol=1),
                          t(as.matrix(vcf[,c("A1", "A2", "C1", "C2"), with=F])))

  gm <- rbind(parent.genomat, f1.genomat)

### make header info
  setnames(vcf, "#CHROM", "chr")

  vcf.ag <- vcf[,list(maxPos=max(POS)), list(chr)]
  setkey(vcf, chr)
  setkey(vcf.ag, chr)

  vcf <- merge(vcf, vcf.ag)
  vcf[,cm:=(POS/maxPos) * 100]


  marker <- matrix(c("marker", vcf$id), nrow=1)
  chr <- matrix(c("chromosome", vcf$chr), nrow=1)
  pos <- matrix(c("pos(cM)", vcf$cm), nrow=1)

  header <- do.call("rbind", list(marker, chr, pos))

  gmh <- rbind(header, gm)


### write output
  writeLines( paste("#founders,",4, sep=""),
               con="/scratch/aob2x/daphnia_hwe_sims/trioPhase/AxC_Rabbit.in"
             )
  options(scipen=999)

   write.table(gmh,
               file="/scratch/aob2x/daphnia_hwe_sims/trioPhase/AxC_Rabbit.in",
               quote=FALSE,
               row.names=FALSE,
               col.names=FALSE,
               sep=",",
               na="N",
               append=TRUE)




  setkey(vcf, "#CHROM")
  founders <- readLines(paste('../01_forward_simulator/', population, '.founders', sep=""))
  vcf <- vcf[.(chromosome), c("#CHROM","POS",founders), with=F]
  # convert to RABBIT format, where ref  = 1, alt=2
  if(any(vcf=="0/0") | any(vcf=="1/1")) {
    # Convert to factor with level 1 = "0/0", level 2 = "1/1"
    vcf[, (founders) := lapply(.SD, factor, levels=c("0/0","1/1")), .SDcols=founders]
    # Convert to numeric
    vcf[, (founders) := lapply(.SD, as.numeric), .SDcols=founders]
    # Subtract 1, such that "0/0" is now 0, "1/1" is now 1, missing is NA
    vcf[, (founders) := lapply(.SD, "-", 1), .SDcols=founders]
    fwrite(vcf, file=paste(population, "/", chromosome, ".RABBIT.vcf" , sep=""))
  }
