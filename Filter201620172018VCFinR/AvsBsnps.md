```
#Load libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(foreach)
        library(SeqArray)
        library(doMC)
        registerDoMC(20)
        library(ggplot2)
        library(tidyverse)

#Load genotype file
        genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
        
#Load sample file        
        load("secondsampstokeep_20190422.Rdata")
        
#Load SNP file
        load("snpsetfilteredformissing_20190423.Rdata")

#Load superclone assignment file
        scs <- fread("Superclones20161718_20190423")

#Pull out A and B individuals 
	AandB <- scs[SC=="A" | SC=="B"]

#Remove males and SM libraries
	AandBnomaleSM <- AandB[clone!="Lab_2019_D8_349Male" & clone!="May_2017_D8_770SM" &
		clone!="Lab_2019_D8_222Male" & clone!="May_2017_D8_773SM"]
	AandBnomaleSMids <- AandBnomaleSM$clone

#Set sequence filters
        seqSetFilter(genofile, sample.id=AandBnomaleSMids)
	seqSetFilter(genofile, variant.id=snpsetfilteredformissing)

#Extract genotypes
	snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		chr = seqGetData(genofile, "chromosome"),
		pos = seqGetData(genofile, "position"),
		dp = seqGetData(genofile, "annotation/info/DP"))

	het <- t(seqGetData(genofile, "$dosage"))
	het <- as.data.table(het)
			
	colnames(het) <- c(seqGetData(genofile, "sample.id"))
	het$variant.ids <- seqGetData(genofile, "variant.id")
			
	setkey(het, variant.ids)
	setkey(snps, variant.ids)
			
	mhet <- merge(snps, het)
			
	sample.ids <- seqGetData(genofile, "sample.id")

#Make long format
	mhetlong <- melt(mhet, measure.vars=sample.ids, variable.name="clone", value.name="dosage")
	setkey(mhetlong, clone)
	setkey(AandBnomaleSM, clone)
	mmhetlong <- merge(AandBnomaleSM, mhetlong)
	
	df2 <- count(mmhetlong, c('variant.ids', 'SC', 'dosage'))
	mmhetlong %>% group_by(variant.ids, SC, dosage) %>% mutate(count = n())

#Make variables for homref, het, and homalt
	mhetlong$homref <- ifelse(mhetlong$dosage=="2", 1, 0)
	mhetlong$het <- ifelse(mhetlong$dosage=="1", 1, 0)
	mhetlong$homalt <- ifelse(mhetlong$dosage=="0", 1, 0)
