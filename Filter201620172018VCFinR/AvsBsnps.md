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

#Transform to wide format
	doscountwide <- dcast(dosagecounts, variant.ids + SC ~ dosage, value.var="N")
	colnames(doscountwide) <- c("variant.ids", "SC", "dos0", "dos1", "dos2")
	yetwider <- dcast(doscountwide, variant.ids~SC, value.var=c("dos0", "dos1", "dos2"))
	max(yetwider$dos0_A, na.rm=TRUE) #88
	max(yetwider$dos0_B, na.rm=TRUE) #29

#Pull out SNPs that are fixed differences between A and B
	A0B2low <- yetwider[dos0_A > 80 & dos2_B > 25] #6381 SNPs
	A0B2mod <- yetwider[dos0_A > 85 & dos2_B > 26] #2515 SNPs
	A0B2high <- yetwider[dos0_A > 86 & dos2_B > 27] #932 SNPs

	A2B0low <- yetwider[dos2_A > 80 & dos0_B > 25] #11261 SNPs
	A2B0mod <- yetwider[dos2_A > 85 & dos0_B > 26] #4360 SNPs
	A2B0high <- yetwider[dos2_A > 86 & dos0_B > 27] #1612 SNPs
	
	ABlow <- rbind(A0B2low, A2B0low) #17642
	ABmod <- rbind(A0B2mod, A2B0mod) #6875
	ABhigh <- rbind(A0B2high, A2B0high) #2544
	
	save(ABlow, file="ABlow_20190426.Rdata")
	save(ABmod, file="ABmod_20190426.Rdata")
	save(ABhigh, file="ABhigh_20190426.Rdata")

```
Let's try looking at these loci in individuals from D8 2018.
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
        

#Load sample file and SNP file       
        load("secondsampstokeep_20190422.Rdata")
	sampsdt <- as.data.table(secondsampstokeep)
	load("ABhigh_20190426.Rdata")
	
#Pull out D8 2018 individuals
	temp <- unlist(strsplit(sampsdt$secondsampstokeep, split="_"))
	mat <- matrix(temp, ncol=4, byrow=TRUE)
	matdat <- as.data.table(mat)
	sampsdt$population <- matdat$V3
	sampsdt$year <- matdat$V2
	sampsdt$season <- matdat$V1
	
	D82018 <- sampsdt$secondsampstokeep[sampsdt$population=="D8" & sampsdt$year=="2018"]
	
	ABhighids <- ABhigh$variant.ids
	
#Set sequence filters
        seqSetFilter(genofile, sample.id=D82018)
	seqSetFilter(genofile, variant.id=ABhighids)

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

 #Count instances of variant.ids/SC/dosage	
	dosagecounts <- mhetlong[, .N, by=list(variant.ids, dosage)]

#Remove NAs
	dosagecounts <- dosagecounts[dosage!="NA"]
	
#Transform to wide format
	doscountwide <- dcast(dosagecounts, variant.ids ~ dosage, value.var="N")
	colnames(doscountwide) <- c("variant.ids", "dos0", "dos1", "dos2")
	doscountwide[is.na(dos0),dos0:=0]
	doscountwide[is.na(dos1),dos1:=0]
	doscountwide[is.na(dos2),dos2:=0]
	doscountwide$total <- doscountwide$dos0+doscountwide$dos1+doscountwide$dos2
	doscountwide$numaltallele <- (doscountwide$dos0*2) + doscountwide$dos1
	doscountwide$propalt <- doscountwide$numaltallele/(doscountwide$total*2)
	doscountwide$foldpropalt <- ifelse(doscountwide$propalt > 0.5, 1-doscountwide$propalt, doscountwide$propalt)
	
	ggplot(doscountwide, aes(x=foldpropalt)) + geom_histogram() + xlim(0,0.5)

#Maybe some but not all are F1 hybrids.
	dosagebyclonedt <- as.data.table(table(mhetlong$clone, mhetlong$dosage))
	colnames(dosagebyclonedt) <- c("clone", "dosage", "N")
	dosagebyclonedtw <- dcast(dosagebyclonedt, clone~dosage, value.var="N")
	
#Need to polarize by A/B
	setkey(ABhigh, variant.ids)
	setkey(mhetlong, variant.ids)
	m <- merge(ABhigh, mhetlong)
	mnoNA <- m[dosage!="NA"]
	mnoNA[is.na(dos0_A),dos0_A:=0]
	mnoNA[is.na(dos1_A),dos1_A:=0]
	mnoNA[is.na(dos2_A),dos2_A:=0]
	mnoNA[is.na(dos0_B),dos0_B:=0]
	mnoNA[is.na(dos1_B),dos1_B:=0]
	mnoNA[is.na(dos2_B),dos2_B:=0]
	mnoNA$dosageAB <- ifelse(mnoNA$dos0_A>0 & mnoNA$dosage==0, 2, ifelse(mnoNA$dos0_A>0 & mnoNA$dosage==2, 0, mnoNA$dosage))
	
#Count instances of variant.ids/SC/dosage	
	dosagecounts <- mnoNA[, .N, by=list(variant.ids, dosageAB)]

#Remove NAs
	dosagecounts <- dosagecounts[dosageAB!="NA"]
	
#Transform to wide format
	doscountwide <- dcast(dosagecounts, variant.ids ~ dosageAB, value.var="N")
	colnames(doscountwide) <- c("variant.ids", "dos0", "dos1", "dos2")
	doscountwide[is.na(dos0),dos0:=0]
	doscountwide[is.na(dos1),dos1:=0]
	doscountwide[is.na(dos2),dos2:=0]
	doscountwide$total <- doscountwide$dos0+doscountwide$dos1+doscountwide$dos2
	doscountwide$numaltallele <- (doscountwide$dos0*2) + doscountwide$dos1
	doscountwide$propalt <- doscountwide$numaltallele/(doscountwide$total*2)
	doscountwide$foldpropalt <- ifelse(doscountwide$propalt > 0.5, 1-doscountwide$propalt, doscountwide$propalt)
	
	ggplot(doscountwide, aes(x=foldpropalt)) + geom_histogram() + xlim(0,0.5)

#Maybe some but not all are F1 hybrids.
	dosagebyclonedt <- as.data.table(table(mnoNA$clone, mnoNA$dosageAB))
	colnames(dosagebyclonedt) <- c("clone", "dosage", "N")
	dosagebyclonedtw <- dcast(dosagebyclonedt, clone~dosage, value.var="N")
	colnames(dosagebyclonedtw) <- c("clone", "dos0", "dos1", "dos2")

#Remove all individuals that are A
	touse <- dosagebyclonedtw[dos2 < 2000]
	setkey(touse, clone)
	setkey(mnoNA, clone)
	mnoNAnoA <- merge(touse, mnoNA)

#Count instances of variant.ids/SC/dosage	
	dosagecounts <- mnoNAnoA[, .N, by=list(variant.ids, dosageAB)]

#Remove NAs
	dosagecounts <- dosagecounts[dosageAB!="NA"]
	
#Transform to wide format
	doscountwide <- dcast(dosagecounts, variant.ids ~ dosageAB, value.var="N")
	colnames(doscountwide) <- c("variant.ids", "dos0", "dos1", "dos2")
	doscountwide[is.na(dos0),dos0:=0]
	doscountwide[is.na(dos1),dos1:=0]
	doscountwide[is.na(dos2),dos2:=0]
	doscountwide$total <- doscountwide$dos0+doscountwide$dos1+doscountwide$dos2
	doscountwide$numaltallele <- (doscountwide$dos0*2) + doscountwide$dos1
	doscountwide$propalt <- doscountwide$numaltallele/(doscountwide$total*2)
	doscountwide$foldpropalt <- ifelse(doscountwide$propalt > 0.5, 1-doscountwide$propalt, doscountwide$propalt)
	
	ggplot(doscountwide, aes(x=foldpropalt)) + geom_histogram() + xlim(0,0.5)

	March2018noAids <- touse$clone
	
	save(March2018noAids, file="March2018noAids_20190426.Rdata")
	
#Try chromosome painting at these SNPs?

#Let's try looking at all individuals from D8 from all years
	
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
        load("ABhigh_20190426.Rdata")
	ABhighids <- ABhigh$variant.ids

#Load superclone assignment file
        scs <- fread("Superclones20161718_20190423")

#Pull out D8 individuals
	temp <- unlist(strsplit(scs$clone, split="_"))
	mat <- matrix(temp, ncol=4, byrow=TRUE)
	matdat <- as.data.table(mat)
	scs$population <- matdat$V3
	scs$year <- matdat$V2
	scs$season <- matdat$V1
	scsD8 <- scs[population=="D8"]

#Remove males and SM libraries
	DnomaleSM <- scsD8[clone!="Lab_2019_D8_349Male" & clone!="May_2017_D8_770SM" &
		clone!="Lab_2019_D8_222Male" & clone!="May_2017_D8_773SM" &
		clone!="May_2017_D8_731SM"]
	DnomaleSMids <- DnomaleSM$clone

#Set sequence filters
        seqSetFilter(genofile, sample.id=DnomaleSMids)
	seqSetFilter(genofile, variant.id=ABhighids)

### get residual heterozygosity
		
	het <- t(seqGetData(genofile, "$dosage"))
	het <- as.data.table(het)
			
	colnames(het) <- c(seqGetData(genofile, "sample.id"))
	het$variant.ids <- seqGetData(genofile, "variant.id")
			
	setkey(het, variant.ids)
	setkey(snps, variant.ids)
			
	mhet <- merge(snps, het)
			
	sample.ids <- seqGetData(genofile, "sample.id")

	mhetlong <- melt(mhet, measure.vars=sample.ids, variable.name="clone", value.name="dosage")

			
			
	normdosagetable <- foreach(i=sample.ids, .combine="rbind")%do%{
		
	data.table(pos=mhet$pos, clone=c(i), normdosage=ifelse(mhet[[i]]==0 & mhet$April_2017_D8_103==0, 2,
		ifelse(mhet[[i]]==2 & mhet$April_2017_D8_103==0, 0, mhet[[i]])))
		
		}


	numclones <- data.table(clone=c(sample.ids), numclone=c(1:176))
	
	setkey(normdosagetable, clone)
	setkey(numclones, clone)
			
	normdosagetablenum <- merge(normdosagetable, numclones)

	numpos <- data.table(pos=c(mhet$pos), numpos=c(1:2544))
	
	setkey(normdosagetablenum, pos)
	setkey(numpos, pos)
			
	normdosagetablenumpos <- merge(normdosagetablenum, numpos)

### Adding in superclones
			
	superclone01 <- fread("Superclones20161718_20190423")

	setkey(normdosagetablenumpos, clone)
	setkey(superclone01, clone)
		
	msc <- merge(normdosagetablenumpos, superclone01)

	msc$numclone <- factor(msc$numclone, levels = msc$numclone[order(msc$SC)])

	msc[, normdosage:=as.factor(normdosage)]
	
	msc$SC <- ifelse(msc$clone=="May_2017_D8_731", "OO", msc$SC)
	
 	temp <- unlist(strsplit(msc$clone, split="_"))
	mat <- matrix(temp, ncol=4, byrow=TRUE)
	matdat <- as.data.table(mat)
	msc$population <- matdat$V3
	msc$year <- matdat$V2
	msc$season <- matdat$V1
	msc$popyearsc <- paste(msc$population,"_",msc$year,"_",msc$SC,sep="")
	msc$numclone <- factor(msc$numclone, levels = msc$numclone[order(msc$popyearsc)])

	ggplot(data=msc, aes(x=numpos, y=numclone, color=normdosage)) + geom_point() + theme(text = element_text(size=7))
