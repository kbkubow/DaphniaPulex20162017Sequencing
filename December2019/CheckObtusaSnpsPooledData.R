#!/usr/bin/env Rscript

### libraries
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(tidyverse)

### First read in SNPs fixed between Obstusa and pulex

    load("/scratch/kbb7sh/Daphnia/MappingDecember2019/snpsPulexObtusafixed_20200121.Rdata")

### Next load in pooled Rdata

    pools <- fread("totalADRDlong")
    colnames(pools) <- c("chr", "pos", "ref", "alt", "RD", "AD", "sample")

### First let's look at D8 where we don't think there should be contamination
    poolsD8 <- pools[sample=="D8Male1" | sample=="D8Male2" | sample=="D8PE1" | sample=="D8PE2"]
    setkey(poolsD8, chr, pos)
    setkey(snpsPulexObtusafixed, chr, pos)
    mD8 <- merge(snpsPulexObtusafixed, poolsD8)

    #See that after merging we retain 4 of 1,283,752 fixed SNPs and 497,993 variable SNPs within the pools

    mD8$propalt <- mD8$AD/(mD8$AD+mD8$RD)
    ggplot(mD8, aes(x=propalt)) + geom_histogram()

    ggplot(data=mD8, aes(x=RD, fill="RD")) + geom_histogram() +
      geom_histogram(aes(x=mD8$AD, fill="AD"))

### How does that compare to DBunk?
    poolsDBunk <- pools[sample=="DBunkMale" | sample=="DBunkPE1" | sample=="DBunkPE2"]
    setkey(poolsDBunk, chr, pos)
    setkey(snpsPulexObtusafixed, chr, pos)
    mDBunk <- merge(snpsPulexObtusafixed, poolsDBunk)

    #See that after merging we retain 4 of 1,283,752 SNPs, the same as with D8

    mDBunk$propalt <- mDBunk$AD/(mDBunk$AD+mDBunk$RD)
    ggplot(mDBunk, aes(x=propalt)) + geom_histogram()


### How does that compare to DCat?
    poolsDCat <- pools[sample=="DCatMale" | sample=="DCatPE1" | sample=="DCatPE2"]
    setkey(poolsDCat, chr, pos)
    setkey(snpsPulexObtusafixed, chr, pos)
    mDCat <- merge(snpsPulexObtusafixed, poolsDCat)

    #See that after merging we retain 4 of 1,283,752 SNPs, the same as with D8

    mDCat$propalt <- mDCat$AD/(mDCat$AD+mDCat$RD)
    ggplot(mDCat, aes(x=propalt)) + geom_histogram()
