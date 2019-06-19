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

### open geno file
		genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### open subset file
    load("subsetindwithcountsctotal_20190501.Rdata")
    
### Pull out pond/year/etc info
    temp <- unlist(strsplit(subsetindwithcountsc$clone, split="_"))
		mat <- matrix(temp, ncol=4, byrow=TRUE)
		matdat <- as.data.table(mat)
		subsetindwithcountsc$population <- matdat$V3
		subsetindwithcountsc$year <- matdat$V2
		subsetindwithcountsc$season <- matdat$V1
    
    DBunk2017April <- subsetindwithcountsc[population=="DBunk" & year=="2017" & season!="May"]
    
    
		
