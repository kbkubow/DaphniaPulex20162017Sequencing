#!/usr/bin/env Rscript

### libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(SeqArray)

### format data properly to work with SNPRelate
        ### subsetted data for PCA stuff
                vcf.fn <- "/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_wOG_90.vcf"
                snpgdsVCF2GDS(vcf.fn, "/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_wOG_90.gds",
                                        method=c("biallelic.only"), snpfirstdim = FALSE)

                seqVCF2GDS(vcf.fn, "/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_wOG_90.seq.gds")


### format data properly to work with SNPRelate
        ### subsetted data for PCA stuff
                vcf.fn <- "/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/GoodSitesB.vcf"
                snpgdsVCF2GDS(vcf.fn, "/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/GoodSitesB.gds",
                                          method=c("biallelic.only"), snpfirstdim = FALSE)

                seqVCF2GDS(vcf.fn, "/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/GoodSitesB.seq.gds")



module load  gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3


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

        genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_wOG_90.seq.gds")

        snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
              chr = seqGetData(genofile, "chromosome"),
              pos = seqGetData(genofile, "position"),
              dp = seqGetData(genofile, "annotation/info/DP"),
              filter = seqGetData(genofile, "annotation/filter"),
              af=seqAlleleFreq(genofile, .progress=T))

        snpsPass <- snps[filter=="PASS"]
        snpsPass99 <- snps[filter=="PASS" | filter=="VQSRTrancheSNP90.00to99.00"]
        snpsPassids <- snpsPass$variant.ids

        genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/GoodSitesB.seq.gds")

        goodsnpsB <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
              chr = seqGetData(genofile, "chromosome"),
              pos = seqGetData(genofile, "position"),
              dp = seqGetData(genofile, "annotation/info/DP"),
              filter = seqGetData(genofile, "annotation/filter"),
              af=seqAlleleFreq(genofile, .progress=T))

### Let's try moving forward with the SNPs that pass the 90 filter

    seqSetFilter(genofile, variant.id=snpsPassids)

    initialsnps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
          chr = seqGetData(genofile, "chromosome"),
          pos = seqGetData(genofile, "position"),
          dp = seqGetData(genofile, "annotation/info/DP"),
          filter = seqGetData(genofile, "annotation/filter"),
          af=seqAlleleFreq(genofile, .progress=T))

    initialsnpsids <- initialsnps$variant.ids


    NsChrRD <- fread("/scratch/kbb7sh/Daphnia/genomefiles/NsandDepthandChrEnd.sorted.500merged.bed")
    colnames(NsChrRD) <- c("chr", "start", "stop")
    NsChrRD$count <- c(1:26531)

    setkey(snps, chr, pos)

    NsChrRDsnps <- foreach(i=NsChrRD$count, .combine="c")%do%{

      c=NsChrRD$chr[[i]]
      s=NsChrRD$start[[i]]
      p=NsChrRD$stop[[i]]

      temp <- snps[J(data.table(chr=c, pos=c(s:p), key="chr,pos")), nomatch=0]
      temp$variant.ids

    }

    save(NsChrRDsnps, file="/scratch/kbb7sh/Daphnia/MappingDecember2019/recalibrate90/wOG/NsChrRDsnps_20200721.Rdata")

    goodsnpsnotinNsChrRD <- setdiff(initialsnpsids, NsChrRDsnps)

    seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRD)

    goodsnpsnotinNsChrRDtable <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
    	chr = seqGetData(genofile, "chromosome"),
    	pos = seqGetData(genofile, "position"),
    	dp = seqGetData(genofile, "annotation/info/DP"))

    save(goodsnpsnotinNsChrRDtable, file="/scratch/kbb7sh/Daphnia/MappingDecember2019/recalibrate90/wOG/goodsnpsnotinNsChrRDtable_woG_20200803.Rdata")

# Started with 1,112,586 SNPs, have 959,705 left after removing NsChrRD SNPs (152,881 SNPs removed)

### Now remove SNPs in repeat masker identified regions
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

    genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_wOG_90.seq.gds")

		load("/scratch/kbb7sh/Daphnia/MappingDecember2019/recalibrate90/wOG/goodsnpsnotinNsChrRDtable_woG_20200803.Rdata")

		RMout <- fread("/scratch/kbb7sh/Daphnia/genomefiles/RMoutHiCGMgoodscaff.bed")
		colnames(RMout) <- c("chr", "start", "stop")
		RMout$count <- c(1:125881)

		setkey(goodsnpsnotinNsChrRDtable, chr, pos)
		initialsnpsB <- goodsnpsnotinNsChrRDtable$variant.ids

		RMoutSNPs <- foreach(i=RMout$count, .combine="c")%do%{

			c=RMout$chr[[i]]
			s=RMout$start[[i]]
			p=RMout$stop[[i]]

			temp <- goodsnpsnotinNsChrRDtable[J(data.table(chr=c, pos=c(s:p), key="chr,pos")), nomatch=0]
			temp$variant.ids

		}

		save(RMoutSNPs, file="/scratch/kbb7sh/Daphnia/MappingDecember2019/recalibrate90/wOG/RMoutSNPs_wOG_20200722.Rdata")

		goodsnpsnotinNsChrRDorRM <- setdiff(initialsnpsB, RMoutSNPs)

		seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRDorRM)

		goodsnpsnotinNsChrRDorRMtable <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
			chr = seqGetData(genofile, "chromosome"),
			pos = seqGetData(genofile, "position"),
			dp = seqGetData(genofile, "annotation/info/DP"))

		save(goodsnpsnotinNsChrRDorRMtable, file="/scratch/kbb7sh/Daphnia/MappingDecember2019/recalibrate90/wOG/goodsnpsnotinNsChrRDorRMtable_wOG_20200722.Rdata")

### Removing repeat masker regions resulted in a loss of 78,966 SNPs with 880,739 SNPs remaining.

###Now remove triallelic snps
			goodsnpsnotinNsChrRDorRM <- goodsnpsnotinNsChrRDorRMtable$variant.ids

			seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRDorRM)

			tri <- (seqGetData(genofile, "$num_allele"))
			tri <- as.data.table(tri)
			tri$variant.ids <- seqGetData(genofile, "variant.id")
			tri$diallelic <- ifelse(tri$tri=="2", 1, 0)

			updatesnpstouse <- tri$variant.ids[tri$tri=="2"]

			save(updatesnpstouse, file="/scratch/kbb7sh/Daphnia/MappingDecember2019/recalibrate90/wOG/updatesnpstouse_depthfiltandnotrialleleic_20200803.Rdata")

			seqSetFilter(genofile, variant.id=updatesnpstouse)

			snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
						chr = seqGetData(genofile, "chromosome"),
						pos = seqGetData(genofile, "position"),
						dp = seqGetData(genofile, "annotation/info/DP"))

### Removing triallelic SNPs resulted in a loss of 37,471 SNPs, with 843,268 remaining

### Overall read depth filtering on snps without removing Repeat Masker

			quantile(snps$dp,probs=c(0.05,0.95))

      #5%   95%
      #3931 12177

			lowRDsnps <- snps[dp < 3931]
			highRDsnps <- snps[dp > 12177]
			lowRDsnps$RD <- c("low")
			highRDsnps$RD <- c("high")
			lowhighRDsnps <- rbind(lowRDsnps, highRDsnps)

			save(lowhighRDsnps, file="/scratch/kbb7sh/Daphnia/MappingDecember2019/recalibrate90/wOG/lowhighRDsnps_20200803.Rdata")

			dpfiltsnps <- snps[dp > 3930 & dp < 12178]
			dpfiltsnpsids <- dpfiltsnps$variant.ids

			save(dpfiltsnps, file="/scratch/kbb7sh/Daphnia/MappingDecember2019/recalibrate90/wOG/dpfiltsnps_20200803.Rdata")

			ggplot(data=dpfiltsnps, aes(x=dp)) + geom_histogram()
			ggplot(data=dpfiltsnps, aes(x=log10(dp))) + geom_histogram()

			seqSetFilter(genofile, variant.id=dpfiltsnpsids)

			snpsG <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
					chr = seqGetData(genofile, "chromosome"),
					pos = seqGetData(genofile, "position"),
					dp = seqGetData(genofile, "annotation/info/DP"))

### Resulted in a removal of 84,283 SNPs, leaving 758,985 SNPs.

### Now let's pull out the already identified pulicaria and obtusa and very low coverage ind, and then do a final check that new individuals are all pulex

		sample.ids <- seqGetData(genofile, "sample.id")
		sampleidsdt <- as.data.table(sample.ids)
		samplestouseB <- sampleidsdt$sample.ids[sample.ids!="Spring_2017_DBunk_340" &
		sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
		sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
		sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
		sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
		sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
		sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
		sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
		sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
		sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43" &
		sample.ids!="2018_Pulicaria_Pond21_22" & sample.ids!="2018_Pulicaria_Pond22_21" &
		sample.ids!="2018_Pulicaria_Pond22_53" & sample.ids!="2018_Pulicaria_Pond22_62" &
		sample.ids!="2018_Pulicaria_Pond22_72"]

		seqSetFilter(genofile, sample.id=samplestouseB)

	### set some global parameters
		maf <- 0.001
		missing.rate <- 0.15
		threads <- 10

	#Let's temporarily use a previous LD SNPset.
		load("./../mhet_20200121.Rdata")
		muhetsub <- data.table(chr=mhet$chr, pos=mhet$pos)
		setkey(snpsG, chr, pos)
		setkey(muhetsub, chr, pos)
		snpstouse <- merge(muhetsub, snpsG)
		snpstouseids <- snpstouse$variant.ids

		pca <- snpgdsPCA(genofile, snp.id=snpstouseids, sample.id=samplestouseB, autosome.only=FALSE, maf=maf,
				missing.rate=missing.rate, num.thread=threads)
		#Working space: 569 samples, 281,919 SNVs

		pc.percent <- pca$varprop*100
		head(round(pc.percent, 2))
		#	[1] 23.16 13.19  8.07  6.78  4.67  4.16

		tab <- data.frame(clone = pca$sample.id,
			EV1 = pca$eigenvect[,1],
			EV2 = pca$eigenvect[,2],
			stringsAsFactors = FALSE)
		head(tab)

		temp <- unlist(strsplit(tab$clone, split="_"))
		mat <- matrix(temp, ncol=4, byrow=TRUE)
		matdat <- as.data.table(mat)
		tab$population <- matdat$V3
		tab$year <- matdat$V2
		tab$season <- matdat$V1
		tab$popseason <- paste(tab$population,"_",tab$season, sep="")
		tab$season2 <- ifelse(tab$season=="April", "Spring", tab$season)
		tab$seasonyear <- paste(tab$year,"_",tab$season2, sep="")

		tab <- as.data.table(tab)

		PCAall <- ggplot(data=tab, aes(x=EV1, y=EV2, color=population)) + geom_point()
		ggsave(PCAall, file="PCAall_20200623.pdf")

		PCAnoW <- ggplot(data=tab[population!="W1" & population!="W6" &
				population!="D10"], aes(x=EV1, y=EV2, color=popseason)) + geom_point()
		ggsave(PCAnoW, file="PCAnoW_20200623.pdf")

### Confirmed! No more non-pulex individuals.

### Ok, now let's remove all SNPs that are fixed within pulex

	sample.ids <- seqGetData(genofile, "sample.id")
	sampleidsdt <- as.data.table(sample.ids)
	samplestouseB <- sampleidsdt$sample.ids[sample.ids!="Spring_2017_DBunk_340" &
	sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
	sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
	sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
	sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
	sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
	sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
	sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
	sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
	sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43" &
	sample.ids!="2018_Pulicaria_Pond21_22" & sample.ids!="2018_Pulicaria_Pond22_21" &
	sample.ids!="2018_Pulicaria_Pond22_53" & sample.ids!="2018_Pulicaria_Pond22_62" &
	sample.ids!="2018_Pulicaria_Pond22_72"]

	seqSetFilter(genofile, sample.id=samplestouseB)

	seqSetFilter(genofile, variant.id=dpfiltsnpsids)

	snps.dt <- data.table(variant.ids=seqGetData(genofile, "variant.id"),
												af=seqAlleleFreq(genofile, .progress=T))

### filter down to ones that are polymorphic in pulex
	pulexPoly <- snps.dt[af!=1 & af!=0]

	setkey(dpfiltsnps, variant.ids)
	setkey(pulexPoly, variant.ids)

	snpsvarPulex <- merge(dpfiltsnps, pulexPoly)

	save(snpsvarPulex, file="/scratch/kbb7sh/Daphnia/MappingDecember2019/recalibrate90/wOG/snpsvarPulex_20200803.Rdata")

### Removing SNPs non variable in D. pulex results in a loss of 398,839 SNPs, with 360,146 remaining.

## Ok, now let's check number of individuals per snp and number of snps per index

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

		genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_90.seq.gds")

		load("snpsvarPulex_20200722.Rdata")
		snpstousevarPulex <- snpsvarPulex$variant.ids
		seqSetFilter(genofile, variant.id=snpstousevarPulex)

		sample.ids <- seqGetData(genofile, "sample.id")
		sampleidsdt <- as.data.table(sample.ids)
		samplestouseB <- sampleidsdt$sample.ids[sample.ids!="Spring_2017_DBunk_340" &
			sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
			sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
			sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
			sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
			sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
			sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
			sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
			sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
			sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43" &
			sample.ids!="2018_Pulicaria_Pond21_22" & sample.ids!="2018_Pulicaria_Pond22_21" &
			sample.ids!="2018_Pulicaria_Pond22_53" & sample.ids!="2018_Pulicaria_Pond22_62" &
			sample.ids!="2018_Pulicaria_Pond22_72"]

	#Pull out genotypes

			seqSetFilter(genofile, sample.id=samplestouseB)

			het <- t(seqGetData(genofile, "$dosage"))
			het <- as.data.table(het)

			colnames(het) <- c(seqGetData(genofile, "sample.id"))
			het$variant.ids <- seqGetData(genofile, "variant.id")

			setkey(het, variant.ids)
			setkey(snpsvarPulex, variant.ids)

			mhet <- merge(snpsvarPulex, het)

	# Set all genotypes to 1

			mhet[mhet == 0] <- 1
			mhet[mhet == 2] <- 1

			mhetlong <- melt(mhet, measure.vars=samplestouseB, variable.name="clone", value.name="dosage")

			save(mhetlong, file="mhetlong_countinggenotypespersnp_20200722.Rdata")

			mhetlong.ag <- mhetlong[,list(numgeno = sum(dosage, na.rm=TRUE)), list(variant.ids) ]

			save(mhetlong.ag, file="mhetlong.ag_countinggenotypespersnp_20200722.Rdata")

  # Having done this, let's look at the distribution of number of individuals genotypes per SNP.

      indpersnp <- ggplot(data=mhetlong.ag, aes(x=numgeno)) + geom_histogram(binwidth=10)
      ggsave(indpersnp, file="indpersnp_20200623.pdf")

    # Let's drop SNPs that are genotyped in less than half the individuals, this drops 1108 SNPs and leaves 439,994
      snpsvarpulexpresentinhalf <- mhetlong.ag$variant.ids[mhetlong.ag$numgeno>282]

      seqSetFilter(genofile, variant.id=snpsvarpulexpresentinhalf)

      save(snpsvarpulexpresentinhalf, file="snpsvarpulexpresentinhalf_20200722.Rdata")

      snpsvarPulex <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
      		chr = seqGetData(genofile, "chromosome"),
      		pos = seqGetData(genofile, "position"),
      		dp = seqGetData(genofile, "annotation/info/DP"))

      write.table(snpsvarPulex, file="snpsvarpulexpresentinhalf_table_20200722", sep="\t", row.names=FALSE, quote=FALSE)

### Let's look for SNPs that are polymorphic within D8, DCat, DBunk

      sample.ids <- seqGetData(genofile, "sample.id")
      sampleidsdt <- as.data.table(sample.ids)

      temp <- unlist(strsplit(as.character(sampleidsdt$sample.ids), split="_"))
      mat <- matrix(temp, ncol=4, byrow=TRUE)
      matdat <- as.data.table(mat)
      sampleidsdt$population <- matdat$V3
      sampleidsdtD8DBunkDCat <- sampleidsdt[population=="D8" | population=="DBunk" | population=="DCat"]

      samplestouseBdt <- sampleidsdtD8DBunkDCat[sample.ids!="Spring_2017_DBunk_340" &
      	sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
      	sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
      	sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
      	sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
      	sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
      	sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
      	sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
      	sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
      	sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43"]

      samplestouseB <- samplestouseBdt$sample.ids

      seqSetFilter(genofile, sample.id=samplestouseB)

      snps.dt <- data.table(variant.ids=seqGetData(genofile, "variant.id"),
      			af=seqAlleleFreq(genofile, .progress=T))

      ### filter down to ones that are polymorphic in pulex and MAF of at least 0.05
      	#pulexD8DCatDBunkPoly <- snps.dt[af!=1 & af!=0]
      	pulexD8DCatDBunkPoly <- snps.dt[af < 0.95 & af>0.05]

      	pulexD8DCatDBunkPolyids <- pulexD8DCatDBunkPoly$variant.ids
      	seqSetFilter(genofile, variant.id=pulexD8DCatDBunkPolyids)

      	save(pulexD8DCatDBunkPolyids, file="pulexD8DCatDBunkPolyids_20200722.Rdata")

      	seqGDS2VCF(genofile, "polyDorsetMAF005recal.vcf.gz")


### Could we also pull out fixed A vs C SNPs (using prior to LD pruned set)
	   sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/Superclones201617182019withObtusaandPulicaria_kingcorr_20200402_wmedrd.txt")
	   scAC <- sc[SC=="A" | SC=="C"]
	   scACids <- scAC$clone
	   seqSetFilter(genofile, sample.id=scACids)

# Pull out genotypes
	het <- t(seqGetData(genofile, "$dosage"))
	het <- as.data.table(het)

	colnames(het) <- c(seqGetData(genofile, "sample.id"))
	het$variant.ids <- seqGetData(genofile, "variant.id")

	setkey(het, variant.ids)
	setkey(snpsvarPulex, variant.ids)

	mhet <- merge(snpsvarPulex, het)

	mhetlong <- melt(mhet, measure.vars=scACids, variable.name="clone", value.name="dosage")

	setkey(scAC, clone)
	setkey(mhetlong, clone)
	mmhetlong <- merge(mhetlong, scAC)
	mmhetlong <- mmhetlong[Nonindependent==0]

	SCcounts <- mmhetlong[, .N, by=list(SC, variant.ids, dosage)]

#Remove NAs
	SCcounts <- SCcounts[dosage!="NA"]
	SCcountsA <- SCcounts[SC=="A"]
	SCcountsC <- SCcounts[SC=="C"]

	SCcountsAwide <- dcast(SCcountsA, variant.ids ~ dosage, value.var="N")
	colnames(SCcountsAwide) <- c("variant.ids", "Ados0", "Ados1", "Ados2")
	SCcountsAwide[is.na(Ados0),Ados0:=0]
	SCcountsAwide[is.na(Ados1),Ados1:=0]
	SCcountsAwide[is.na(Ados2),Ados2:=0]

	SCcountsCwide <- dcast(SCcountsC, variant.ids ~ dosage, value.var="N")
	colnames(SCcountsCwide) <- c("variant.ids", "Cdos0", "Cdos1", "Cdos2")
	SCcountsCwide[is.na(Cdos0),Cdos0:=0]
	SCcountsCwide[is.na(Cdos1),Cdos1:=0]
	SCcountsCwide[is.na(Cdos2),Cdos2:=0]

	setkey(SCcountsAwide, variant.ids)
	setkey(SCcountsCwide, variant.ids)
	SCcountsAC <- merge(SCcountsAwide, SCcountsCwide)

	SCcountsAC$Atot <- SCcountsAC$Ados0+SCcountsAC$Ados1+SCcountsAC$Ados2
	SCcountsAC$Ctot <- SCcountsAC$Cdos0+SCcountsAC$Cdos1+SCcountsAC$Cdos2

	SCcountsACfilt <- SCcountsAC[Atot>59 & Ctot>19]

	SCcountsACfiltACfix <- SCcountsACfilt[Ados0==Atot & Cdos2==Ctot]
	SCcountsACfiltCAfix <- SCcountsACfilt[Ados2==Atot & Cdos0==Ctot]

	ACfixed <- rbind(SCcountsACfiltACfix, SCcountsACfiltCAfix)
	ACfixedids <- ACfixed$variant.ids
	save(ACfixed, file="ACfixed_20200722.Rdata")

	Ahet <- SCcountsACfilt[Ados1==Atot]
	save(Ahet, file="Ahet_20200722.Rdata")
	Ahetids <- Ahet$variant.ids

	seqSetFilter(genofile, variant.id=Ahetids)
	seqSetFilter(genofile, sample.id=c("AxB_R1_P110_A", "D8179_R1_P58_A"))

	het <- t(seqGetData(genofile, "$dosage"))
	het <- as.data.table(het)

	colnames(het) <- c(seqGetData(genofile, "sample.id"))
	het$variant.ids <- seqGetData(genofile, "variant.id")

	setkey(het, variant.ids)
	setkey(snpsvarPulex, variant.ids)

	mhet <- merge(snpsvarPulex, het)

	Chet <- SCcountsACfilt[Cdos1==Ctot]
	save(Chet, file="Chet_20200722.Rdata")
	Chetids <- Chet$variant.ids

	scB <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
	selfedC <- scB[SC=="selfedC"]
	selfedCids <- selfedC$clone

	seqSetFilter(genofile, variant.id=Chetids)
	seqSetFilter(genofile, sample.id=selfedCids)

	het <- t(seqGetData(genofile, "$dosage"))
	het <- as.data.table(het)

	colnames(het) <- c(seqGetData(genofile, "sample.id"))
	het$variant.ids <- seqGetData(genofile, "variant.id")

	setkey(het, variant.ids)
	setkey(snpsvarPulex, variant.ids)

	mhet <- merge(snpsvarPulex, het)

	mhetlong <- melt(mhet, measure.vars=selfedCids, variable.name="clone", value.name="dosage")

	SCcounts <- mhetlong[, .N, by=list(clone, dosage)]

	SCcounts <- SCcounts[dosage!="NA"]

	SCcountswide <- dcast(SCcounts, clone ~ dosage, value.var="N")
	colnames(SCcountswide) <- c("variant.ids", "dos0", "dos1", "dos2")
	SCcountswide[is.na(dos0),dos0:=0]
	SCcountswide[is.na(dos1),dos1:=0]
	SCcountswide[is.na(dos2),dos2:=0]

	SCcountswide$tot <- SCcountswide$dos0+SCcountswide$dos1+SCcountswide$dos2
	SCcountswide$prophet <- SCcountswide$dos1/SCcountswide$tot


#Now let's look at these SNPs in everyone

	seqResetFilter(genofile)

	sample.ids <- seqGetData(genofile, "sample.id")
	sampleidsdt <- as.data.table(sample.ids)
	samplestouseB <- sampleidsdt$sample.ids[sample.ids!="Spring_2017_DBunk_340" &
		sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
		sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
		sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
		sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
		sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
		sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
		sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
		sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
		sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43" &
		sample.ids!="2018_Pulicaria_Pond21_22" & sample.ids!="2018_Pulicaria_Pond22_21" &
		sample.ids!="2018_Pulicaria_Pond22_53" & sample.ids!="2018_Pulicaria_Pond22_62" &
		sample.ids!="2018_Pulicaria_Pond22_72"]


	seqSetFilter(genofile, sample.id=samplestouseB)
	seqSetFilter(genofile, variant.id=ACfixedids)



# Get genotypes

	het <- t(seqGetData(genofile, "$dosage"))
	het <- as.data.table(het)

	colnames(het) <- c(seqGetData(genofile, "sample.id"))
	het$variant.ids <- seqGetData(genofile, "variant.id")

	setkey(het, variant.ids)
	setkey(snpsvarPulex, variant.ids)

	mhet <- merge(snpsvarPulex, het)

	mhetlong <- melt(mhet, measure.vars=samplestouseB, variable.name="clone", value.name="dosage")

	setkey(mhetlong, variant.ids)
	setkey(ACfixed, variant.ids)
	mhetlongACinfo <- merge(mhetlong, ACfixed)

	mhetlongACinfo$ACsnp <- ifelse(mhetlongACinfo$Ados0 > mhetlongACinfo$Cdos2, "Asnp", "Csnp")
	mhetlongACinfo$dosageB <- ifelse(mhetlongACinfo$dosage==1, 1, ifelse(mhetlongACinfo$dosage==0 & mhetlongACinfo$ACsnp=="Asnp", 2, ifelse(
		mhetlongACinfo$dosage==2 & mhetlongACinfo$ACsnp=="Asnp", 0, mhetlongACinfo$dosage
		)))

	clonecounts <- mhetlongACinfo[, .N, by=list(clone, dosageB)]
	clonecounts <- clonecounts[dosageB!="NA"]

	clonecountswide <- dcast(clonecounts, clone ~ dosageB, value.var="N")
	colnames(clonecountswide) <- c("clone", "HomC", "Het", "HomA")
	clonecountswide[is.na(HomC),HomC:=0]
	clonecountswide[is.na(Het),Het:=0]
	clonecountswide[is.na(HomA),HomA:=0]

	clonecountswide$total <- clonecountswide$HomC+clonecountswide$Het+clonecountswide$HomA
	clonecountswide$propC <- clonecountswide$HomC/clonecountswide$total
	clonecountswide$propHet <- clonecountswide$Het/clonecountswide$total
	clonecountswide$propA <- clonecountswide$HomA/clonecountswide$total

	setkey(clonecountswide, clone)
	setkey(sc, clone)
	mclonecountswide <- merge(clonecountswide, sc, all.x=TRUE)

	mclonecountswide$ACF1hybrid <- ifelse(mclonecountswide$propHet > 0.9, 1, 0)

	save(mclonecountswide, file="mclonecountswide_ACfixedSNPs_20200722.Rdata")


  ### Let's make a vcf for KING with all D8, DBunk, and DCat individuals, and the lab generated clones

  		sample.ids <- seqGetData(genofile, "sample.id")
  		sampleidsdt <- as.data.table(sample.ids)

  		temp <- unlist(strsplit(as.character(sampleidsdt$sample.ids), split="_"))
  		mat <- matrix(temp, ncol=4, byrow=TRUE)
  		matdat <- as.data.table(mat)
  		sampleidsdt$population <- matdat$V3
  		sampleidsdt$population <- ifelse(sampleidsdt$population=="Dcat", "DCat", sampleidsdt$population)
  		sampleidsdtsub <- sampleidsdt[population!="D10" & population!="DLily" & population!="DMud" &
  			population!="DOil" & population!="Dramp" & population!="W1" & population!="W6"]

  		samplestouseBdt <- sampleidsdtsub[sample.ids!="Spring_2017_DBunk_340" &
  		sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
  			sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
  			sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
  			sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
  			sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
  			sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
  			sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
  			sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
  			sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43" &
  			sample.ids!="2018_Pulicaria_Pond21_22" & sample.ids!="2018_Pulicaria_Pond22_21" &
  			sample.ids!="2018_Pulicaria_Pond22_53" & sample.ids!="2018_Pulicaria_Pond22_62" &
  			sample.ids!="2018_Pulicaria_Pond22_72" & sample.ids!="April_2017_D8_515R" &
  			sample.ids!="Lab_2019_D8_222Male" & sample.ids!="Lab_2019_D8_349Male" &
  			sample.ids!="May_2017_D8_731SM" & sample.ids!="May_2017_D8_770SM" &
  			sample.ids!="May_2017_D8_773SM" & sample.ids!="Spring_2017_DBunk_116SM" &
  			sample.ids!="Spring_2017_DBunk_347SM" & sample.ids!="Spring_2017_DBunk_73SM" &
  			sample.ids!="Spring_2016_D8_8.1" & sample.ids!="March_2018_D8_18030" &
  			sample.ids!="March_2018_DCat_18004"]

  		samplestouseB <- samplestouseBdt$sample.ids

  		seqSetFilter(genofile, sample.id=samplestouseB)

  		seqSetFilter(genofile, variant.id=pulexD8DCatDBunkPolyids)

  		seqGDS2VCF(genofile, "polyDorsetMAF005withlabgeneratedrecal.vcf.gz")

#For now lets proceed with SNP pruning and IBS superclone assignment, and see where low coverage individuals fall out
# For now setting maf pretty low, so SNPs are kept even if present in 1 individual, can also try in 2 ind to compare

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
        library(dplyr)
        library(tidyverse)

#Load genotype file
  genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_90.seq.gds")

#Load SNP file
  load("snpsvarpulexpresentinhalf_20200722.Rdata")
  seqSetFilter(genofile, variant.id=snpsvarpulexpresentinhalf)

# Set individual filter
  sample.ids <- seqGetData(genofile, "sample.id")
  sampleidsdt <- as.data.table(sample.ids)

  temp <- unlist(strsplit(as.character(sampleidsdt$sample.ids), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  sampleidsdt$population <- matdat$V3
  sampleidsdt$population <- ifelse(sampleidsdt$population=="Dcat", "DCat", sampleidsdt$population)
  sampleidsdtsub <- sampleidsdt[population=="D8" | population=="DBunk" | population=="DCat" |
    population=="D10" | population=="DLily" | population=="DMud" | population=="DOil" |
    population=="Dramp" | population=="W1" | population=="W6"]

  samplestouseBdt <- sampleidsdtsub[sample.ids!="Spring_2017_DBunk_340" &
  sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
  sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
  sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
  sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
  sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
  sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
  sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
  sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
  sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43" &
  sample.ids!="2018_Pulicaria_Pond21_22" & sample.ids!="2018_Pulicaria_Pond22_21" &
  sample.ids!="2018_Pulicaria_Pond22_53" & sample.ids!="2018_Pulicaria_Pond22_62" &
  sample.ids!="2018_Pulicaria_Pond22_72" & sample.ids!="April_2017_D8_515R" &
  sample.ids!="Lab_2019_D8_222Male" & sample.ids!="Lab_2019_D8_349Male" &
  sample.ids!="May_2017_D8_731SM" & sample.ids!="May_2017_D8_770SM" &
  sample.ids!="May_2017_D8_773SM" & sample.ids!="Spring_2017_DBunk_116SM" &
  sample.ids!="Spring_2017_DBunk_347SM" & sample.ids!="Spring_2017_DBunk_73SM" &
  sample.ids!="Spring_2016_D8_8.1" & sample.ids!="March_2018_D8_18030" &
  sample.ids!="March_2018_DCat_18004"]

  samplestouseB <- samplestouseBdt$sample.ids

  seqSetFilter(genofile, sample.id=samplestouseB)

### set some global parameters
  maf <- 0.001
  missing.rate <- 0.15
  threads <- 10

### SNP pruning (removing lab generated clones)
  set.seed(10000)
  snpset01 <- snpgdsLDpruning(genofile, snp.id=snpsvarpulexpresentinhalf, sample.id=samplestouseB,
  autosome.only=FALSE, maf=maf,missing.rate=missing.rate, slide.max.bp=500, ld.threshold=0.1)
  finalsetsnpset01 <-unlist(snpset01[c(1:62)])
  finalsetsnpset01dt <- as.data.table(finalsetsnpset01)

  save(finalsetsnpset01, file="finalsetsnpset01_20200722.Rdata")
  #This resulted in 112,511 SNPs

  seqSetFilter(genofile, variant.id=finalsetsnpset01)

  snpsvarPulexLDprune <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
      chr = seqGetData(genofile, "chromosome"),
      pos = seqGetData(genofile, "position"),
      dp = seqGetData(genofile, "annotation/info/DP"))

  write.table(snpsvarPulexLDprune, file="finalsetsnpset01pulex_table_20200722", sep="\t", row.names=FALSE, quote=FALSE)

  ### IBS

  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  setkey(sc, clone)
  sctouse <- sc[Nonindependent==0 & LabGenerated==0]
  sctouseids <- sctouse$clone

  #First lets try IBS with we keep SNPs present in 1 individuals
  		### set some global parameters
  			maf <- 0.001
  			missing.rate <- 0.15
  			threads <- 10

  		ibs <- snpgdsIBS(genofile, snp.id=finalsetsnpset01, sample.id=sctouseids, num.thread=20, maf=maf,
  					missing.rate=0.15, autosome.only = FALSE)

  	### a bit of re-formating of the ibs matrix
  		ibs.mat <- ibs$ibs
  		rownames(ibs.mat) <- ibs$sample.id
  		colnames(ibs.mat) <- ibs$sample.id

  	### make the IBs matrix long form
  		ibs.matdt <- as.data.table(ibs.mat)
  		ibs.matdt$cloneA <- sctouseids
  		ibs.long<- melt(ibs.matdt, measure.vars=sctouseids, variable.name="cloneB", value.name="IBS")

  		ibs.long <- na.omit(ibs.long)

  	### Identify cutoff by graphing

  		ggplot(data=ibs.long, aes(x=IBS)) + geom_histogram(binwidth=0.001)

  		ggplot(data=ibs.long[IBS > 0.92 & IBS < 0.99], aes(x=IBS)) + geom_histogram(binwidth=0.001)

  		constrainedIBSdistB <- ggplot(data=ibs.long[IBS > 0.95 & IBS < 0.97], aes(x=IBS)) + geom_histogram(binwidth=0.001)

      ggplot(data=ibs.long[IBS > 0.92 & IBS < 0.98 & cloneA!="April_2017_D8_214" & cloneB!="April_2017_D8_214" &
        cloneA!="Spring_2016_D8_8.3" & cloneB!="Spring_2016_D8_8.3" &
        cloneA!="Spring_2016_D8_8.12" & cloneB!="Spring_2016_D8_8.12"], aes(x=IBS)) + geom_histogram(binwidth=0.001)

  		ggsave(constrainedIBSdistB, file="constrainedIBSdistB_20200623.pdf")

  		### Identify superclones
  			### give temporary labels to super clones based on identity of clone B
  					superclone <- ibs.long[IBS>.94]
  					superclone[,SC.sub := as.numeric(as.factor(cloneB))]

  			### collapse nested superclones
  					superclone.o <- foreach(sc.i=unique(superclone$SC.sub), .combine="c")%do%{
  							paste(sort(unique(superclone[cloneB%in%superclone[SC.sub==sc.i]$cloneA]$cloneA)), collapse=";")
  							}

  					superclone.o.uniq <- unique(superclone.o)

  					sc.dt <- foreach(i=superclone.o.uniq, .combine="rbind")%do%{
  								data.table(clone=strsplit(i, ";")[[1]],
  								superClone.size=length(strsplit(i, ";")[[1]]),
  								superClone.index=which(i==superclone.o.uniq))
  					}
  					sc.dt[,superClone.sizeRank := as.numeric(as.factor(rank(-superClone.size, ties="average")))]

  					sc.dt <- as.data.table(sc.dt %>%
  							mutate(SCnum = group_indices_(sc.dt, .dots=c("superClone.sizeRank", "superClone.index"))))

  					#temp <- unlist(strsplit(sc.dt$clone, split="_"))
  					#mat <- matrix(temp, ncol=4, byrow=TRUE)
  					#matdat <- as.data.table(mat)
  					#sc.dt$population <- matdat$V3
  					#sc.dt$year <- matdat$V2

  					save(sc.dt, file="tmp_sc.dt_20200722.Rdata")


### Let's try looking at resequenced clones
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
        library(dplyr)
        library(tidyverse)

#Load genotype file
  genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_90.seq.gds")

#Load SNP file
  load("snpsvarpulexpresentinhalf_20200722.Rdata")
  seqSetFilter(genofile, variant.id=snpsvarpulexpresentinhalf)

# Load SC file

  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  scC <- sc[SC=="C" & Nonindependent==0]
  scCids <- scC$clone
  seqSetFilter(genofile, sample.id=scCids)

  Csnps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
        chr = seqGetData(genofile, "chromosome"),
        pos = seqGetData(genofile, "position"),
        dp = seqGetData(genofile, "annotation/info/DP"),
        af=seqAlleleFreq(genofile, .progress=T))

  # Get genotypes

  	het <- t(seqGetData(genofile, "$dosage"))
  	het <- as.data.table(het)

  	colnames(het) <- c(seqGetData(genofile, "sample.id"))
  	het$variant.ids <- seqGetData(genofile, "variant.id")

  	setkey(het, variant.ids)
  	setkey(Csnps, variant.ids)

  	mhet <- merge(Csnps, het)

  	mhetlong <- melt(mhet, measure.vars=scCids, variable.name="clone", value.name="dosage")

    variantcounts <- mhetlong[, .N, by=list(variant.ids, dosage)]
    variantcounts <- variantcounts[dosage!="NA"]

    variantcountswide <- dcast(variantcounts, variant.ids ~ dosage, value.var="N")
    colnames(variantcountswide) <- c("variant.id", "dos0", "dos1", "dos2")
    variantcountswide[is.na(dos0),dos0:=0]
    variantcountswide[is.na(dos1),dos1:=0]
    variantcountswide[is.na(dos2),dos2:=0]

    variantcountswide$total <- variantcountswide$dos0+variantcountswide$dos1+variantcountswide$dos2
    variantcountswide$propalt <- ((variantcountswide$dos0*2)+variantcountswide$dos1)/(variantcountswide$total*2)
    variantcountswide$oneover2n <- ifelse(variantcountswide$dos1==1 & variantcountswide$dos2==variantcountswide$total-1 |
      variantcountswide$dos1==1 & variantcountswide$dos0==variantcountswide$total-1, 1, 0)
    variantcountswide$foldedaf <- ifelse(variantcountswide$propalt > 0.5, 1-variantcountswide$propalt, variantcountswide$propalt)

    Coneover2nsnps <- variantcountswide[oneover2n==1]
    Coneover2nsnpsids <- Coneover2nsnps$variant.id

    seqSetFilter(genofile, variant.id=Coneover2nsnpsids)
    seqSetFilter(genofile, sample.id=c("April_2017_D8_515R", "May_2017_D8_515"))

    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(Csnps, variant.ids)

    mhet <- merge(Csnps, het)

    seqSetFilter(genofile, variant.id=Coneover2nsnpsids)
    seqSetFilter(genofile, sample.id=c("Lab_2019_D8_222Male", "April_2017_D8_222"))

    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(Csnps, variant.ids)

    mhet <- merge(Csnps, het)


    scA <- sc[SC=="A" & Nonindependent==0]
    scAids <- scA$clone
    seqSetFilter(genofile, sample.id=scAids)

    Asnps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
          chr = seqGetData(genofile, "chromosome"),
          pos = seqGetData(genofile, "position"),
          dp = seqGetData(genofile, "annotation/info/DP"),
          af=seqAlleleFreq(genofile, .progress=T))

    # Get genotypes

    	het <- t(seqGetData(genofile, "$dosage"))
    	het <- as.data.table(het)

    	colnames(het) <- c(seqGetData(genofile, "sample.id"))
    	het$variant.ids <- seqGetData(genofile, "variant.id")

    	setkey(het, variant.ids)
    	setkey(Asnps, variant.ids)

    	mhet <- merge(Asnps, het)

    	mhetlong <- melt(mhet, measure.vars=scAids, variable.name="clone", value.name="dosage")

      variantcounts <- mhetlong[, .N, by=list(variant.ids, dosage)]
      variantcounts <- variantcounts[dosage!="NA"]

      variantcountswide <- dcast(variantcounts, variant.ids ~ dosage, value.var="N")
      colnames(variantcountswide) <- c("variant.id", "dos0", "dos1", "dos2")
      variantcountswide[is.na(dos0),dos0:=0]
      variantcountswide[is.na(dos1),dos1:=0]
      variantcountswide[is.na(dos2),dos2:=0]

      variantcountswide$total <- variantcountswide$dos0+variantcountswide$dos1+variantcountswide$dos2
      variantcountswide$propalt <- ((variantcountswide$dos0*2)+variantcountswide$dos1)/(variantcountswide$total*2)
      variantcountswide$oneover2n <- ifelse(variantcountswide$dos1==1 & variantcountswide$dos2==variantcountswide$total-1 |
        variantcountswide$dos1==1 & variantcountswide$dos0==variantcountswide$total-1, 1, 0)

      Aoneover2nsnps <- variantcountswide[oneover2n==1]
      Aoneover2nsnpsids <- Aoneover2nsnps$variant.id

      seqSetFilter(genofile, variant.id=Aoneover2nsnpsids)
      seqSetFilter(genofile, sample.id=c("Lab_2019_D8_349Male", "April_2017_D8_349"))

      het <- t(seqGetData(genofile, "$dosage"))
      het <- as.data.table(het)

      colnames(het) <- c(seqGetData(genofile, "sample.id"))
      het$variant.ids <- seqGetData(genofile, "variant.id")

      setkey(het, variant.ids)
      setkey(Csnps, variant.ids)

      mhet <- merge(Csnps, het)
