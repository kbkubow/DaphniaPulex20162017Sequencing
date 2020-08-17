#module load intel/18.0  intelmpi/18.0 gsl/2.4 R/3.6.3 boost/1.60.0; R
 # module load

### libraries
	library(data.table)
	library(foreach)
	library(SeqArray)
  library(LEA)

setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/")

### open GDS file
	genofile <- seqOpen("MapJune2020_ann.seq.gds", readonly=TRUE)

### snpFilter file
  snpFilter <- fread("snpsvarpulexpresentinhalf_table_20200623")

### load metadata file
  samps <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  samps[,set:=paste(Species, population, sep="_")]

  samps.ag.sc <- samps[!SC%in%c("OO", "AxCF1", "selfedA", "selfedC") & Nonindependent==0,
                    list(clone=clone[which.max(medrd)]), SC]
  samps.ag.oo <- samps[SC%in%c("OO"), list(clone=clone), SC]

  samps.ag <- rbind(samps.ag.sc, samps.ag.oo)

### write vcf file as input for LEA
  seqSetFilter(genofile, variant.id=snpFilter$variant.ids, sample.id=samps.ag$clone)

  seqGDS2VCF(genofile, "/scratch/aob2x/daphnia_hwe_sims/daps4lea.vcf", info.var=character(0), fmt.var=character(0), use_Rsamtools=TRUE,
    verbose=TRUE)

### lea convert to GENO object
  vcf2geno("/scratch/aob2x/daphnia_hwe_sims/daps4lea.vcf",
            output.file = "/scratch/aob2x/daphnia_hwe_sims/daps4lea.geno",
            force = TRUE)

### ancestry test
snmf("/scratch/aob2x/daphnia_hwe_sims/daps4lea.geno",
    K=c(1:10),
    project = "continue",
    repetitions = 1, CPU = 1,
    alpha = 10, tolerance = 0.00001, entropy = FALSE, percentage = 0.05, I=10000, iterations = 200, ploidy = 2, seed = -1)
