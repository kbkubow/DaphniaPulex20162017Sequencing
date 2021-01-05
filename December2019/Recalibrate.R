#!/usr/bin/env Rscript

### libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(SeqArray)

### format data properly to work with SNPRelate
        ### subsetted data for PCA stuff
                vcf.fn <- "/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_90.vcf"
                snpgdsVCF2GDS(vcf.fn, "/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_90.gds",
                                        method=c("biallelic.only"), snpfirstdim = FALSE)

                seqVCF2GDS(vcf.fn, "/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_90.seq.gds")


### format data properly to work with SNPRelate
        ### subsetted data for PCA stuff
                vcf.fn <- "/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_99.vcf"
                snpgdsVCF2GDS(vcf.fn, "/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_99.gds",
                                          method=c("biallelic.only"), snpfirstdim = FALSE)

                seqVCF2GDS(vcf.fn, "/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_99.seq.gds")

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

        genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_90.seq.gds")

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

    save(NsChrRDsnps, file="NsChrRDsnps_20200721.Rdata")

    goodsnpsnotinNsChrRD <- setdiff(initialsnpsids, NsChrRDsnps)

    seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRD)

    goodsnpsnotinNsChrRDtable <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
    	chr = seqGetData(genofile, "chromosome"),
    	pos = seqGetData(genofile, "position"),
    	dp = seqGetData(genofile, "annotation/info/DP"))

    save(goodsnpsnotinNsChrRDtable, file="goodsnpsnotinNsChrRDtable_20200722.Rdata")

# Started with 1,423,120 SNPs, have 1,224,102 left after removing NsChrRD SNPs (199,018 SNPs removed)

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

		genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann_90.seq.gds")

		load("goodsnpsnotinNsChrRDtable_20200722.Rdata")

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

		save(RMoutSNPs, file="RMoutSNPs_20200722.Rdata")

		goodsnpsnotinNsChrRDorRM <- setdiff(initialsnpsB, RMoutSNPs)

		seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRDorRM)

		goodsnpsnotinNsChrRDorRMtable <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
			chr = seqGetData(genofile, "chromosome"),
			pos = seqGetData(genofile, "position"),
			dp = seqGetData(genofile, "annotation/info/DP"))

		save(goodsnpsnotinNsChrRDorRMtable, file="goodsnpsnotinNsChrRDorRMtable_20200722.Rdata")

### Removing repeat masker regions resulted in a loss of 100,659 SNPs with 1,123,443 SNPs remaining.

###Now remove triallelic snps
			goodsnpsnotinNsChrRDorRM <- goodsnpsnotinNsChrRDorRMtable$variant.ids

			seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRDorRM)

			tri <- (seqGetData(genofile, "$num_allele"))
			tri <- as.data.table(tri)
			tri$variant.ids <- seqGetData(genofile, "variant.id")
			tri$diallelic <- ifelse(tri$tri=="2", 1, 0)

			updatesnpstouse <- tri$variant.ids[tri$tri=="2"]

			save(updatesnpstouse, file="updatesnpstouse_depthfiltandnotrialleleic_20200722.Rdata")

			seqSetFilter(genofile, variant.id=updatesnpstouse)

			snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
						chr = seqGetData(genofile, "chromosome"),
						pos = seqGetData(genofile, "position"),
						dp = seqGetData(genofile, "annotation/info/DP"))

### Removing triallelic SNPs resulted in a loss of 44,974 SNPs, with 1,078,469 remaining

### Overall read depth filtering on snps without removing Repeat Masker

			quantile(snps$dp,probs=c(0.05,0.95))

      #5%   95%
      #3827 12311

			lowRDsnps <- snps[dp < 3827]
			highRDsnps <- snps[dp > 12311]
			lowRDsnps$RD <- c("low")
			highRDsnps$RD <- c("high")
			lowhighRDsnps <- rbind(lowRDsnps, highRDsnps)

			save(lowhighRDsnps, file="lowhighRDsnps_20200722.Rdata")

			dpfiltsnps <- snps[dp > 3826 & dp < 12312]
			dpfiltsnpsids <- dpfiltsnps$variant.ids

			save(dpfiltsnps, file="dpfiltsnps_20200722.Rdata")

			ggplot(data=dpfiltsnps, aes(x=dp)) + geom_histogram()
			ggplot(data=dpfiltsnps, aes(x=log10(dp))) + geom_histogram()

			seqSetFilter(genofile, variant.id=dpfiltsnpsids)

			snpsG <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
					chr = seqGetData(genofile, "chromosome"),
					pos = seqGetData(genofile, "position"),
					dp = seqGetData(genofile, "annotation/info/DP"))

### Resulted in a removal of 107,780 SNPs, leaving 970,689 SNPs.

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

	save(snpsvarPulex, file="snpsvarPulex_20200722.Rdata")

### Removing SNPs non variable in D. pulex results in a loss of 529,587 SNPs, with 441,102 remaining.

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
