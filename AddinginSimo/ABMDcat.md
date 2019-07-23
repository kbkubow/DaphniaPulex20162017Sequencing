#Load libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(foreach)
        library(SeqArray)
        library(doMC)
        registerDoMC(20)
        library(ggplot2)
        library(viridis)
        library(tidyverse)

# Load genotype file and snp and sample filter files

        genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

        load("finalsnpstousewSimoids_20190430.Rdata")
        
        load("secondsampstokeepwSimo_20190430.Rdata")
        
# Set sequence filters

        seqSetFilter(genofile, sample.id=secondsampstokeepwSimo)
        seqSetFilter(genofile, variant.id=finalsnpstousewSimoids)
        
# Load superclone file
      
        scs <- fread("Superclones20161718updated_20190501")

# Pull out all A, B, M, amd DCat18004 clones

        scsABM <- scs[superClone.index=="2" | superClone.index=="3" | superClone.index=="33" | superClone.index=="86"]

#Remove males and SM libraries

        scsAnomaleSM <- scsABM[clone!="Lab_2019_D8_349Male" & clone!="May_2017_D8_770SM"]
        scsAnomaleSMids <- scsAnomaleSM$clone
        
# Reset sequence filter

        seqSetFilter(genofile, sample.id=scsAnomaleSMids)
