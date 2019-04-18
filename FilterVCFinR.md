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
                vcf.fn <- "/mnt/spicy_3/Karen/Sp2017/NewMapping/totalnewmapBfiltsnps10bpindels_snps_filter_pass_lowGQmiss.vcf"
                snpgdsVCF2GDS(vcf.fn, "/mnt/spicy_3/Karen/Sp2017/NewMapping/totalnewmapBfiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds", 
                                        method=c("biallelic.only"), snpfirstdim = FALSE)
        
                seqVCF2GDS(vcf.fn, "/mnt/spicy_3/Karen/Sp2017/NewMapping/totalnewmapBfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
'''
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

		genofile <- seqOpen("/mnt/spicy_3/Karen/Sp2017/NewMapping//totalnewmapBfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

		snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
					chr = seqGetData(genofile, "chromosome"),
					pos = seqGetData(genofile, "position"),
					dp = seqGetData(genofile, "annotation/info/DP"))

    	ggplot(data=snps, aes(x=dp)) + geom_histogram()
	ggplot(data=snps, aes(x=log10(dp))) + geom_histogram()
```
![SNPsInitialReadDepthDistribution](/https://github.com/kbkubow/DaphniaPulex20162017Sequencing/blob/master/SNPsInitialReadDepthDistribution.tiff)

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

	save(NsChrRDsnps, file="NsChrRDsnps_20190418.Rdata")

	seqSetFilter(genofile, variant.id=NsChrRDsnps) # 331,095
		
	NsChrRDsnpssnps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		chr = seqGetData(genofile, "chromosome"),
		pos = seqGetData(genofile, "position"),
		dp = seqGetData(genofile, "annotation/info/DP"))

	ggplot(data=NsChrRDsnpssnps, aes(x=dp)) + geom_histogram()		
	ggplot(data=NsChrRDsnpssnps, aes(x=log10(dp))) + geom_histogram()		
		
	goodsnpsnotinNsChrRD <- setdiff(initialsnps, NsChrRDsnps)
