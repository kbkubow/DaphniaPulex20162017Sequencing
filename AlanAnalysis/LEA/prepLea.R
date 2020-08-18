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
#snmf("/scratch/aob2x/daphnia_hwe_sims/daps4lea.geno",
#    K=c(1:10),
#    project = "continue",
#    repetitions = 1, CPU = 1,
#    alpha = 10, tolerance = 0.00001, entropy = FALSE, percentage = 0.05, I=10000, iterations = 200, ploidy = 2, seed = -1)

#obj.snmf = snmf("/scratch/aob2x/daphnia_hwe_sims/daps4lea.geno", K = 3, alpha = 100, project = "new")



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

#scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/Q_samps.Rdata ~/Q_samps.Rdata
# R
library(LEA)
library(foreach)
library(data.table)
library(ggplot2)
library(pophelper)
library(cowplot); theme_set(theme_cowplot())
load("~/Q_samps.Rdata")

processQ <- function(k, samp) {
	qi<-as.data.table(q.list[[k]])
	qi <- as.data.table(qi)
	qi[,gr:=apply(qi, 1, which.max)]
	qi[,samp:=samp]
	qil <- melt(qi, measure.vars=paste("V", 1:k, sep=""))

	qil.ag <- foreach(gr.i=unique(qil$gr), .combine="rbind")%do%{
		qil.sub <- qil[gr==gr.i][variable==paste("V", gr.i, sep="")]
		qil.sub[,ord:=rank(value, ties.method="first")]
		qil.sub
	}

	setkey(qil, samp)
	setkey(qil.ag, samp)

	qil <- merge(qil, qil.ag)
	qil[,ordf:=factor(ord, levels=sort(unique(ord)))]
	setnames(qil, "samp", "clone")
	qil
}


k <- processQ(k=5, samp=samps$clone)
k <- merge(k, samps, by="clone")


ggplot(k, aes(x=ord, y=value.x, fill=as.factor(variable.x))) +
geom_bar(position="stack", stat="identity") +
facet_grid(~gr.y, scale="free_x", space="free") +
geom_text(data=k[Species=="pulicaria"], aes(x=ord, y=1, label="*"), size=2) +
geom_text(data=k[Species=="obtusa"], aes(x=ord, y=1, label="$"), size=2)
