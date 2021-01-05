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


  scJ <- sc[SC=="J" & clone!="March20_2018_DBunk_5"]
  scJids <- scJ$clone

  seqSetFilter(genofile, sample.id=scJids)

  # Pull out genotypes
  	het <- t(seqGetData(genofile, "$dosage"))
  	het <- as.data.table(het)

  	colnames(het) <- c(seqGetData(genofile, "sample.id"))
  	het$variant.ids <- seqGetData(genofile, "variant.id")

  	setkey(het, variant.ids)
  	setkey(snpsvarPulex, variant.ids)

  	mhet <- merge(snpsvarPulex, het)

  	mhetlong <- melt(mhet, measure.vars=scJids, variable.name="clone", value.name="dosage")

  	variantcounts <- mhetlong[, .N, by=list(variant.ids, dosage)]

  #Remove NAs
  	variantcounts <- variantcounts[dosage!="NA"]

  	variantcountswide <- dcast(variantcounts, variant.ids ~ dosage, value.var="N")
  	colnames(variantcountswide) <- c("variant.ids", "dos0", "dos1", "dos2")
  	variantcountswide[is.na(dos0),dos0:=0]
  	variantcountswide[is.na(dos1),dos1:=0]
  	variantcountswide[is.na(dos2),dos2:=0]
    variantcountswide$total <- variantcountswide$dos0+variantcountswide$dos1+variantcountswide$dos2
    Jhet <- variantcountswide[dos1==total]
    Jhetids <- Jhet$variant.ids

    seqSetFilter(genofile, variant.id=Jhetids)

    seqSetFilter(genofile, sample.id=c("March20_2018_DBunk_5"))

    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snpsvarPulex, variant.ids)

    mhet <- merge(snpsvarPulex, het)
