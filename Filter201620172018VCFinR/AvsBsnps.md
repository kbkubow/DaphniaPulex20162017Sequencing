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
	
#Count instances of variant.ids/SC/dosage	
	dosagecounts <- mmhetlong[, .N, by=list(variant.ids, SC, dosage)]

#Remove NAs
	dosagecounts <- dosagecounts[dosage!="NA"]
	save(dosagecounts, file="dosagecounts_AB_20190426.Rdata")
	doscountwide <- dcast(dosagecounts, variant.ids + SC ~ dosage, value.var="N")
	colnames(doscountwide) <- c("variant.ids", "SC", "dos0", "dos1", "dos2")
	yetwider <- dcast(doscountwide, variant.ids~SC, value.var=c("dos0", "dos1", "dos2"))
	max(yetwider$dos0_A, na.rm=TRUE) #88
	max(yetwider$dos0_B, na.rm=TRUE) #29
	
	A0B2low <- yetwider[dos0_A > 80 & dos2_B > 25] #6381 SNPs
	A0B2mod <- yetwider[dos0_A > 85 & dos2_B > 26] #2515 SNPs
	A0B2high <- yetwider[dos0_A > 86 & dos2_B > 27] #932 SNPs

	A2B0low <- yetwider[dos2_A > 80 & dos0_B > 25] #11261 SNPs
	A2B0mod <- yetwider[dos2_A > 85 & dos0_B > 26] #4360 SNPs
	A2B0high <- yetwider[dos2_A > 86 & dos0_B > 27] #1612 SNPs
	
	ABlow <- rbind(A0B2low, A2B0low) #17642
	ABmod <- rbind(A0B2mod, A2B0mod) #6875
	ABhigh <- rbind(A0B2high, A2B0high) #2544
	
	

	
