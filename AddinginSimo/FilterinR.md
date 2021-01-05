Now I will import the VCF into R, and do additional filtering there.
First, I will convert the VCF into gds and seq.gds formats.
```
#!/usr/bin/env Rscript

### libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(SeqArray)
        library(tidyverse)

### format data properly to work with SNPRelate
        ### subsetted data for PCA stuff
                vcf.fn <- "/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.vcf"
                snpgdsVCF2GDS(vcf.fn, "/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds",
                                        method=c("biallelic.only"), snpfirstdim = FALSE)

                seqVCF2GDS(vcf.fn, "/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
```
Now let's open the file and take an intial loook at the SNPs.
```
### libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(SeqArray)
        library(tidyverse)

		genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

		snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
					chr = seqGetData(genofile, "chromosome"),
					pos = seqGetData(genofile, "position"),
					dp = seqGetData(genofile, "annotation/info/DP"))

    ggplot(data=snps, aes(x=dp)) + geom_histogram()
	  ggplot(data=snps, aes(x=log10(dp))) + geom_histogram()
```
Now let's filter out SNPs that occur in areas flagged as having too high or low read depth when mapping the D84A 10X Illumina short reads to the D84Aimages reference genome. Will also filter out SNPs on the edges of runs of Ns and at the ends of scaffolds.
```
	NsChrRD <- fread("/mnt/spicy_3/Karen/RefGenome/Dovetail/HiCnew/NsandDepthandChrEnd.sorted.500merged.bed")
	colnames(NsChrRD) <- c("chr", "start", "stop")
	NsChrRD$count <- c(1:26531)

	setkey(snps, chr, pos)
	initialsnps <- snps$variant.ids

	NsChrRDsnps <- foreach(i=NsChrRD$count, .combine="c")%do%{

	c=NsChrRD$chr[[i]]
	s=NsChrRD$start[[i]]
	p=NsChrRD$stop[[i]]

	temp <- snps[J(data.table(chr=c, pos=c(s:p), key="chr,pos")), nomatch=0]
	temp$variant.ids

	}

	save(NsChrRDsnps, file="NsChrRDsnpswMarch2018wSimo_20190430.Rdata")

	seqSetFilter(genofile, variant.id=NsChrRDsnps)

	NsChrRDsnpssnps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		chr = seqGetData(genofile, "chromosome"),
		pos = seqGetData(genofile, "position"),
		dp = seqGetData(genofile, "annotation/info/DP"))

	ggplot(data=NsChrRDsnpssnps, aes(x=dp)) + geom_histogram()		
	ggplot(data=NsChrRDsnpssnps, aes(x=log10(dp))) + geom_histogram()		

	goodsnpsnotinNsChrRD <- setdiff(initialsnps, NsChrRDsnps)

	seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRD)

	goodsnpsnotinNsChrRDtable <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		chr = seqGetData(genofile, "chromosome"),
		pos = seqGetData(genofile, "position"),
		dp = seqGetData(genofile, "annotation/info/DP"))

	ggplot(data=goodsnpsnotinNsChrRDtable, aes(x=dp)) + geom_histogram()
	ggplot(data=goodsnpsnotinNsChrRDtable, aes(x=log10(dp))) + geom_histogram()

	save(goodsnpsnotinNsChrRDtable, file="goodsnpsnotinNsChrRDtablewMarch2018wSimo_20190430.Rdata")


```
481,148 SNPs were removed based on this filtering. 2,404,522 SNPs remain.

Now lets filter out snps based on overall read depth. Removing SNPs that are in the tails (too low and too high) of the distribution.
```
### Overall read depth filtering on snps without removing Repeat Masker

	quantile(goodsnpsnotinNsChrRDtable$dp,probs=c(0.015,0.975))

 #1.5% 97.5%
 #2034 11692

	lowRDsnps <- goodsnpsnotinNsChrRDtable[dp < 2034]
	highRDsnps <- goodsnpsnotinNsChrRDtable[dp > 11692]
	lowRDsnps$RD <- c("low")
	highRDsnps$RD <- c("high")
	lowhighRDsnps <- rbind(lowRDsnps, highRDsnps)

	save(lowhighRDsnps, file="lowhighRDsnps_20200102.Rdata")

	dpfiltsnps <- goodsnpsnotinNsChrRDtable[dp > 2033 & dp < 11693]
	dpfiltsnpsids <- dpfiltsnps$variant.ids

	save(dpfiltsnps, file="dpfiltsnps_20200102.Rdata")

	ggplot(data=dpfiltsnps, aes(x=dp)) + geom_histogram()
	ggplot(data=dpfiltsnps, aes(x=log10(dp))) + geom_histogram()

```
This depth filtering resulted in a removal of 96,157 SNPs, such that 2,308,365 SNPs remain. Here is the distribution of read depth of the remaining SNPs.

Now let's filter SNPs that have high or low read depth for specific superclones.

```
### Superclone specific depth filtering

	seqSetFilter(genofile, variant.id=dpfiltsnpsids)

	superclones <- fread("/mnt/spicy_3/Karen/Sp2017/supercloneassignment20180816.csv")
	colnames(superclones) <- c("clone", "pond", "year", "season", "SC094", "Plate", "SCnew", "SCnewer", "notes")
	superclones$clone <- ifelse(superclones$clone=="May_2017_DBunk_549", "May_2017_DBunk_549B", superclones$clone)

	Aclones <- superclones$clone[superclones$SC094=="A"]
	Bclones <- superclones$clone[superclones$SC094=="B"]
	Dclones <- superclones$clone[superclones$SC094=="D"]
	Hclones <- superclones$clone[superclones$SC094=="H"]
	Iclones <- superclones$clone[superclones$SC094=="I"]
	Lclones <- superclones$clone[superclones$SC094=="L"]
	Mclones <- superclones$clone[superclones$SC094=="M"]

	supercloneinput <- c("A", "B", "D", "H", "I", "L", "M")

	computingoutofdpsnps <- foreach(x=1:length(supercloneinput), .combine="rbind")%do%{

	tempgroup <- supercloneinput[x]
	tempgroupids <- superclones$clone[superclones$SC094==tempgroup]

	seqSetFilter(genofile, sample.id=tempgroupids)

### per-strain
    ### calculate depth

        dp <- t(seqGetData(genofile, "annotation/format/DP")$data)
        dp <- as.data.table(dp)
        dp[,variant.ids:=dpfiltsnpsids]

		sampleids <- as.data.table(seqGetData(genofile, "sample.id"))
		variantidhead <- as.data.table("variant.ids")
		totalids <- rbind(sampleids, variantidhead)
		colnames(dp) <- c(totalids$V1)

	    setkey(dp, variant.ids)
	    setkey(dpfiltsnps, variant.ids)

	    m <- merge(dp, dpfiltsnps)

		readdepthrow <- melt(m, measure.vars=sampleids, variable.name="clone", value.name="dp")

		readdepthrow.ag <- readdepthrow[,list(totrd=as.double(sum(dp1, na.rm=TRUE))),
					 list(variant.ids)]

		setkey(readdepthrow.ag, variant.ids)
		setkey(dpfiltsnps, variant.ids)
		mreaddepthrow.ag <- merge(readdepthrow.ag, dpfiltsnps)

		ggplot(data=mreaddepthrow.ag, aes(x=totrd)) + geom_histogram()

		tempquant <- as.data.table(quantile(mreaddepthrow.ag$totrd,probs=c(0.025,0.975)))


  		mreaddepthrow.aghigh <- mreaddepthrow.ag[totrd > tempquant$V1[[2]]]
  		mreaddepthrow.aghigh$highlow <- c("high")
   		mreaddepthrow.aglow <- mreaddepthrow.ag[totrd < tempquant$V1[[1]]]
   		mreaddepthrow.aglow$highlow <- c("low")
		mreaddepthrow.aglownottouse <- rbind(mreaddepthrow.aghigh, mreaddepthrow.aglow)
		mreaddepthrow.aglownottouse$clone <- tempgroup
		mreaddepthrow.aglownottouse$numind <- length(tempgroupids)
		mreaddepthrow.aglownottouse

	}

	save(computingoutofdpsnps, file="computingoutofdpsnps_wMarch2018wSImo_20190430.Rdata")

	computingoutofdpsnpsids <- data.table(variant.ids=computingoutofdpsnps$variant.ids)
	computingoutofdpsnpsidsuniq <- unique(computingoutofdpsnpsids)

	totalnotusescdepth <- computingoutofdpsnpsidsuniq$variant.ids

	dpfiltandscfilt <- setdiff(dpfiltsnpsids, totalnotusescdepth)

	seqSetFilter(genofile, variant.id=dpfiltandscfilt)

	snpsdpfiltscfilt <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		chr = seqGetData(genofile, "chromosome"),
		pos = seqGetData(genofile, "position"),
		dp = seqGetData(genofile, "annotation/info/DP"))

	save(dpfiltandscfilt, file="dpfiltandscfilt_wMarch2018wSimo_20190430.Rdata")

	ggplot(data=snpsdpfiltscfilt, aes(x=dp)) + geom_histogram()
```
This superclone specific depth filtering resulted in 824,421 SNPs being removed, leaving 2,057,306 SNPs. Here is the read depth distribution of the remaining SNPs.

Now lets remove trialleleic SNPs.

```
      	library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(SeqArray)
        library(tidyverse)

	genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

	load(file="dpfiltandscfilt_wMarch2018wSimo_20190430.Rdata")

	seqSetFilter(genofile, variant.id=dpfiltandscfilt)

	tri <- (seqGetData(genofile, "$num_allele"))		
	tri <- as.data.table(tri)
	tri$variant.ids <- seqGetData(genofile, "variant.id")
	tri$diallelic <- ifelse(tri$tri=="2", 1, 0)

	updatesnpstouse <- tri$variant.ids[tri$tri=="2"]

	save(updatesnpstouse, file="updatesnpstouse_depthfiltandnotrialleleic_20200102.Rdata")

	seqSetFilter(genofile, variant.id=updatesnpstouse)

#	seqGDS2VCF(genofile, "updatesnpstouse_depthfiltandnotrialleleic_20190419.vcf")

```
Filtering for triallelic SNPs removed 35,640 SNPs. 2,368,882 SNPs remain.

Now let's try to figure out what Samples and SNPs should be dropped to minimize missing data.
```
	seqSetFilter(genofile, variant.id=updatesnpstouse)

	snpsG <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		chr = seqGetData(genofile, "chromosome"),
		pos = seqGetData(genofile, "position"),
		dp = seqGetData(genofile, "annotation/info/DP"))

#Pull out genotypes

	het <- t(seqGetData(genofile, "$dosage"))
	het <- as.data.table(het)

	colnames(het) <- c(seqGetData(genofile, "sample.id"))
	het$variant.ids <- seqGetData(genofile, "variant.id")

	setkey(het, variant.ids)
	setkey(snpsG, variant.ids)

	mhet <- merge(snpsG, het)

# Set all genotypes to 1

	mhet[mhet == 0] <- 1
	mhet[mhet == 2] <- 1

	save(mhet, file="mhetgenosas1wSimo_20190430.Rdata")

# Sum genotypes across individuals, to get number of genotype sites per individual

	genocounts <- as.data.table(colSums(Filter(is.numeric, mhet), na.rm=TRUE))

	genocounts <- genocounts[4:435]

	sample.ids <- seqGetData(genofile, "sample.id")

	genocounts$sample.ids <- sample.ids

	colnames(genocounts) <- c("numsitesgenotype", "sample.ids")

	setkey(genocounts, numsitesgenotype)

	save(genocounts, file="genocountswSimo_20190430.Rdata")

# Let's look at the distribution of number of genotype sites across individuals. Haven't gone beyond this point. Both Simo are really low.

	ggplot(data=genocounts, aes(x=numsitesgenotype)) + geom_histogram(binwidth=10000)

	sampsbelow1450000<- genocounts$sample.ids[genocounts$numsitesgenotype<1450000]
	sampsabove1450000 <- genocounts$sample.ids[genocounts$numsitesgenotype>1450000]

	save(sampsbelow1450000, file="sampsbelow1450000wSimo_20190430.Rdata")
	save(sampsabove1450000, file="sampsabove1450000wSimo_20190430.Rdata")

```
For now I am dropping all individuals with sites genotyped below 1,375,000. This results in dropping 17 individuals. March20_2018_DBunk_39, March20_2018_D8_19, Spring_2017_D8Simo_1, Fall_2016_D10_54, Spring_2017_DRampsSimo_3, April_2017_DBunk_23, Spring_2017_DBunk_19, Spring_2016_D8_8.10, March20_2018_DBunk_21, March20_2018_D8_8, March20_2018_D8_3, Spring_2017_DBunk_27, March20_2018_D8_16, March20_2018_DBunk_1, March20_2018_DBunk_22, March20_2018_DBunk_37, and Spring_2017_DBunk_76.

```
# Removing all individuals with sites genotyped below 1,375,000 from the mhet file.

	mhet[, (sampsbelow1450000) := NULL]

```
Now lets look at the distribution of number of individuals genotyped per SNP.

```
# Having done this, let's look at the distribution of number of individuals genotypes per SNP.

	mhetlong <- melt(mhet, measure.vars=sampsabove1450000, variable.name="clone", value.name="dosage")

	mhetlong.ag <- mhetlong[,list(numgeno = sum(dosage, na.rm=TRUE)), list(variant.ids) ]

	#save(mhetlong.ag, file="mhetlong.agwSimo_20190430.Rdata")

	ggplot(data=mhetlong.ag, aes(x=numgeno)) + geom_histogram(binwidth=10)

```

For now I will drop SNPs that are genotyped in less than 350 individuals (353/415 = 85% individuals genotyped). This results in dropping 204,357 SNPs, with 1,804,302 SNPs remaining.
```
	mhetlong.ag353snpsids <- mhetlong.ag$variant.ids[mhetlong.ag$numgeno>353]

	mhetlong.ag353snpsidsdt <- data.table(variant.ids=mhetlong.ag353snpsids)

	finalsnpstousewSimo <- mhetlong.ag353snpsidsdt
	finalsnpstousewSimoids <- mhetlong.ag353snpsids

	save(finalsnpstousewSimo, file="finalsnpstousewSimo_20190430.Rdata")
	save(finalsnpstousewSimoids, file="finalsnpstousewSimoids_20190430.Rdata")

#	mhet[, (sampsbelow1375000) := NULL]

	setkey(mhetlong.ag353snpsids, variant.ids)
	setkey(mhet, variant.ids)
	mmhet <- merge(mhetlong.ag353snpsids, mhet)

	seqSetFilter(genofile, sample.id=sampsabove1450000)

	genocountsB <- as.data.table(colSums(Filter(is.numeric, mmhet), na.rm=TRUE))

	genocountsB <- genocountsB[4:418]

	sample.ids <- seqGetData(genofile, "sample.id")

	genocountsB$sample.ids <- sample.ids

	colnames(genocountsB) <- c("numsitesgenotype", "sample.ids")

	setkey(genocountsB, numsitesgenotype)

	ggplot(data=genocountsB, aes(x=numsitesgenotype)) + geom_histogram(binwidth=10000)

```
![DistofNumSitesGenoPerIndAfter1stDrop](DistofNumSitesGenoPerIndAfter1stDrop.tiff)
I am going to drop individuals that have fewer than 1,443,442 (80% of total) SNPs genotyped. This results in dropping 5 additional individuals: Spring_2016_D8_8.20, Spring_2017_DBunk_73SM, March20_2018_D8_5, March20_2018_D8_29, and May_2017_DBunk_509. This leaves 410 individuals.

```
	secondsampstodrop <- genocountsB$sample.ids[genocountsB$numsitesgenotype<1443442]
	secondsampstokeep <- genocountsB$sample.ids[genocountsB$numsitesgenotype>1443442]

	save(secondsampstodrop, file="secondsampstodropwSimo_20190430.Rdata")
	save(secondsampstokeep, file="secondsampstokeepwSimo_20190430.Rdata")

	mmhet[, (secondsampstodrop) := NULL]

```
DID NOT GET PAST HERE
Now let's LD prune these SNPs, then run IBS and do superclone identification. We know D. obtusa are in this data set. Let's leave them in for now, identify them in the IBS superclone analysis, and then we can remove them and see if that changes things.
So first LD pruning.
```
  #Load libraries
  	library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(SeqArray)
        library(tidyverse)

# Load genotyping and filtered snp file
	genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

	load("finalsnpstousewSimoids_20190430.Rdata")
	load("secondsampstokeepwSimo_20190430.Rdata")

	seqSetFilter(genofile, variant.id=finalsnpstousewSimoids)
	seqSetFilter(genofile, sample.id=secondsampstokeep)


### set some global parameters
		maf <- 0.01
		missing.rate <- 0.15
		threads <- 10

### SNP pruning
		set.seed(10000)
		snpset01 <- snpgdsLDpruning(genofile, snp.id=finalsnpstousewSimoids, sample.id=secondsampstokeep,
		autosome.only=FALSE, maf=maf,missing.rate=missing.rate, slide.max.bp=500, ld.threshold=0.1)

		finalsetsnpset01 <-unlist(snpset01[c(1:62)])
		finalsetsnpset01dt <- as.data.table(finalsetsnpset01)

		save(finalsetsnpset01, file="finalsetsnpset01wSimo_20190430.Rdata")
		save(finalsetsnpset01dt, file="finalsetsnpset01dtwSimo_20190430.Rdata")

```
LD pruning results in retention of 154,731 SNPs.
