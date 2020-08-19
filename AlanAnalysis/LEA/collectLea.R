#module load intel/18.0  intelmpi/18.0 gsl/2.4 R/3.6.3 boost/1.60.0; R

#### post process
library(LEA)
library(foreach)
library(data.table)

load(file="/scratch/aob2x/daphnia_hwe_sims/dap.snmf.Rdata")




q.list <- foreach(k=1:10)%do%{
	Q(dap.snmf, K = k)
}


setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/")

### load metadata file
  samps <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  samps[,set:=paste(Species, population, sep="_")]




	samps.vcf <- fread("/scratch/aob2x/daphnia_hwe_sims/daps4lea.vcf", skip="#CHROM", nrows=1)
	samps.vcf <- names(samps.vcf)[-(1:9)]
	setkey(samps, clone)

	samps <- samps[J(samps.vcf)]
	samps[,x:=1:dim(samps)[1]]

save(q.list, samps, file="/scratch/aob2x/Q_samps.Rdata")
