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
![SNPsInitialReadDepthDistribution](SNPsInitialReadDepthDistribution.tiff)
