#!/usr/bin/env Rscript

    ### libraries
      library(data.table)
      library(foreach)
      library(ggplot2)
      library(tidyverse)
      library(doMC)
      registerDoMC(20)

### load files

  load("totalpropsimACother_20200818.Rdata")
  load("totalpropsimDBunkD8DCatvsOut_20200818.Rdata")

### Manipulate files
  totalgenomewide <- rbind(totalpropsimACother,totalpropsimDBunkD8DCatvsOut)

  totalgenomewidenoog <- totalgenomewide[comppond!="AvsC_A_obtusa" & comppond!="AvsC_A_pulicaria" &
    comppond!="AvsC_C_obtusa" & comppond!="AvsC_C_pulicaria" & comppond!="AvsC_D8DBunkDCatvsObtusa" &
    comppond!="AvsC_D8DBunkDCatvsPulicaria"]

  totalgenomewidenoognoAC <- totalgenomewidenoog[comppond!="AvsC_A_D10" & comppond!="AvsC_A_W1" &
    comppond!="AvsC_A_W6" & comppond!="AvsC_C_D10" & comppond!="AvsC_C_W1" & comppond!="AvsC_C_W6"]

  totalgenomewidenoognoAC$type <- ifelse(totalgenomewidenoognoAC$comppond=="AvsC_D8D8", "AvsC_D8D8", ifelse(
    totalgenomewidenoognoAC$comppond=="AvsC_DBunkDBunk" | totalgenomewidenoognoAC$comppond=="AvsC_DCatDCat",
    "AvsC_WithinPond", ifelse(totalgenomewidenoognoAC$comppond=="AvsC_D8DBunk" |
    totalgenomewidenoognoAC$comppond=="AvsC_D8DCat" | totalgenomewidenoognoAC$comppond=="AvsC_DBunkDCat",
    "AvsC_BetweenPond", totalgenomewidenoognoAC$comppond
  )))

  totalgenomewidenoognoAC$type <- factor(totalgenomewidenoognoAC$type, levels=c("AvsC_D8D8", "AvsC_WithinPond",
    "AvsC_BetweenPond", "AvsC_D8DBunkDCatvsD10", "AvsC_D8DBunkDCatvsW1", "AvsC_D8DBunkDCatvsW6"))

  totalgenomewidenoognoACsub <- totalgenomewidenoognoAC[, c("comppond", "mean_ratiosim", "sd_ratiosim", "type"), with=FALSE]

  save(totalgenomewidenoognoACsub, file="totalgenomewidenoognoACsub.Rdata")
