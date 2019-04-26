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
        AandB <- scs[SC=="A" | SC=="B"]

#Set sequence filters
        seqSetFilter(genofile, sample.id=secondsampstokeep)
	seqSetFilter(genofile, variant.id=snpsetfilteredformissing)
