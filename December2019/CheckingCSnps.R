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


  load("pulexD8DCatDBunkPolyids_20200722.Rdata")

  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  seqSetFilter(genofile, variant.id=pulexD8DCatDBunkPolyids)

  snpsvarPulex <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
      chr = seqGetData(genofile, "chromosome"),
      pos = seqGetData(genofile, "position"),
      dp = seqGetData(genofile, "annotation/info/DP"))


  scC <- sc[SC=="C" & Nonindependent==0]
  scCids <- scC$clone

  seqSetFilter(genofile, sample.id=scCids)

  # Pull out genotypes
  	het <- t(seqGetData(genofile, "$dosage"))
  	het <- as.data.table(het)

  	colnames(het) <- c(seqGetData(genofile, "sample.id"))
  	het$variant.ids <- seqGetData(genofile, "variant.id")

  	setkey(het, variant.ids)
  	setkey(snpsvarPulex, variant.ids)

  	mhet <- merge(snpsvarPulex, het)

  	mhetlong <- melt(mhet, measure.vars=scCids, variable.name="clone", value.name="dosage")

  	variantcounts <- mhetlong[, .N, by=list(variant.ids, dosage)]

  #Remove NAs
  	variantcounts <- variantcounts[dosage!="NA"]

  	variantcountswide <- dcast(variantcounts, variant.ids ~ dosage, value.var="N")
  	colnames(variantcountswide) <- c("variant.ids", "dos0", "dos1", "dos2")
  	variantcountswide[is.na(dos0),dos0:=0]
  	variantcountswide[is.na(dos1),dos1:=0]
  	variantcountswide[is.na(dos2),dos2:=0]
    variantcountswide$total <- variantcountswide$dos0+variantcountswide$dos1+variantcountswide$dos2
    variantcountswide$altalleles <- variantcountswide$dos1 + (2*variantcountswide$dos0)
    C1over2N <- variantcountswide[altalleles==(total*2)-1]
    C2over2N <- variantcountswide[altalleles==(total*2)-2]
    C1over2Nids <- C1over2N$variant.ids
    C2over2Nids <- C2over2N$variant.ids


    seqSetFilter(genofile, variant.id=C1over2Nids)
    #9,831 SNPs

    seqSetFilter(genofile, sample.id=c("April_2017_D8_222", "May_2017_D8_515"))

    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snpsvarPulex, variant.ids)

    mhet <- merge(snpsvarPulex, het)

    D8222_1over2N <- mhet[April_2017_D8_222==1]
    D8515_1over2N <- mhet[May_2017_D8_515==1]

    seqSetFilter(genofile, variant.id=C2over2Nids)
    #2,035 SNPs

    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snpsvarPulex, variant.ids)

    mhet <- merge(snpsvarPulex, het)

    D8222_2over2N <- mhet[April_2017_D8_222==1]
    D8515_2over2N <- mhet[May_2017_D8_515==1]

    seqSetFilter(genofile, variant.id=D8222_2over2N$variant.ids)

    selfedCs <- sc[SC=="selfedC"]
    selfedCsids <- selfedCs$clone

    seqSetFilter(genofile, sample.id=(c(selfedCsids, "Lab_2019_D8_222Male")))

    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snpsvarPulex, variant.ids)

    mhet <- merge(snpsvarPulex, het)

    mhetlong <- melt(mhet, measure.vars=c(selfedCsids, "Lab_2019_D8_222Male"), variable.name="clone", value.name="dosage")

    variantcounts <- mhetlong[, .N, by=list(clone, dosage)]

  #Remove NAs
    variantcounts <- variantcounts[dosage!="NA"]

    variantcountswide <- dcast(variantcounts, clone~ dosage, value.var="N")
    colnames(variantcountswide) <- c("variant.ids", "dos0", "dos1")
    variantcountswide[is.na(dos0),dos0:=0]
    variantcountswide[is.na(dos1),dos1:=0]
    variantcountswide$total <- variantcountswide$dos0+variantcountswide$dos1

    seqSetFilter(genofile, variant.id=D8515_2over2N$variant.ids)

    selfedCs <- sc[SC=="selfedC"]
    selfedCsids <- selfedCs$clone

    seqSetFilter(genofile, sample.id=(c(selfedCsids, "April_2017_D8_515R")))

    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snpsvarPulex, variant.ids)

    mhet <- merge(snpsvarPulex, het)

    mhetlong <- melt(mhet, measure.vars=c(selfedCsids, "April_2017_D8_515R"), variable.name="clone", value.name="dosage")

    variantcounts <- mhetlong[, .N, by=list(clone, dosage)]

  #Remove NAs
    variantcounts <- variantcounts[dosage!="NA"]

    variantcountswide <- dcast(variantcounts, clone~ dosage, value.var="N")
    colnames(variantcountswide) <- c("variant.ids", "dos0", "dos1")
    variantcountswide[is.na(dos0),dos0:=0]
    variantcountswide[is.na(dos1),dos1:=0]
    variantcountswide$total <- variantcountswide$dos0+variantcountswide$dos1


    seqSetFilter(genofile, variant.id=D8222_1over2N$variant.ids)

    selfedCs <- sc[SC=="selfedC"]
    selfedCsids <- selfedCs$clone

    seqSetFilter(genofile, sample.id=(c(selfedCsids, "Lab_2019_D8_222Male")))

    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snpsvarPulex, variant.ids)

    mhet <- merge(snpsvarPulex, het)

    mhetlong <- melt(mhet, measure.vars=c(selfedCsids, "Lab_2019_D8_222Male"), variable.name="clone", value.name="dosage")

    variantcounts <- mhetlong[, .N, by=list(clone, dosage)]

  #Remove NAs
    variantcounts <- variantcounts[dosage!="NA"]

    variantcountswide <- dcast(variantcounts, clone~ dosage, value.var="N")
    colnames(variantcountswide) <- c("variant.ids", "dos0", "dos1")
    variantcountswide[is.na(dos0),dos0:=0]
    variantcountswide[is.na(dos1),dos1:=0]
    variantcountswide$total <- variantcountswide$dos0+variantcountswide$dos1

    seqSetFilter(genofile, variant.id=D8515_1over2N$variant.ids)

    selfedCs <- sc[SC=="selfedC"]
    selfedCsids <- selfedCs$clone

    seqSetFilter(genofile, sample.id=(c(selfedCsids, "April_2017_D8_515R")))

    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snpsvarPulex, variant.ids)

    mhet <- merge(snpsvarPulex, het)

    mhetlong <- melt(mhet, measure.vars=c(selfedCsids, "April_2017_D8_515R"), variable.name="clone", value.name="dosage")

    variantcounts <- mhetlong[, .N, by=list(clone, dosage)]

  #Remove NAs
    variantcounts <- variantcounts[dosage!="NA"]

    variantcountswide <- dcast(variantcounts, clone~ dosage, value.var="N")
    colnames(variantcountswide) <- c("variant.ids", "dos0", "dos1")
    variantcountswide[is.na(dos0),dos0:=0]
    variantcountswide[is.na(dos1),dos1:=0]
    variantcountswide$total <- variantcountswide$dos0+variantcountswide$dos1
