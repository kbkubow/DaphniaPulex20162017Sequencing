#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(SeqArray)

### load parents and format parental haplotypes
  vcf <- fread("rQTL_quartet.consensus.header.noNA.vcf", skip="#CHROM", na.strings="./.")
  vcf[,A1:=as.numeric(tstrsplit(A, "\\||\\/")[[1]]) + 1]
  vcf[,A2:=as.numeric(tstrsplit(A, "\\||\\/")[[2]]) + 1]
  vcf[,C1:=as.numeric(tstrsplit(C, "\\||\\/")[[1]]) + 1]
  vcf[,C2:=as.numeric(tstrsplit(C, "\\||\\/")[[2]]) + 1]

### load AxC F1s
  f1s <- fread("scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.delim")


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
