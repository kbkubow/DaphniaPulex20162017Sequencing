#!/usr/bin/env Rscript

### libraries
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(tidyverse)
        library(lme4)


### Load data file

  meso <- fread("SampleCountsMeso2019Week7.csv")
  meso$CloneRep <- paste(meso$Clone, meso$Replicate, sep="_")
  meso$SC <- ifelse(meso$Clone=="D8179" | meso$Clone=="D8349", "A", ifelse(
    meso$Clone=="D8222" | meso$Clone=="D8515", "C", "AxC"
  ))

### Cleaning up and adding new variables
  meso$Clone <- as.factor(meso$Clone)
  meso[is.na(JuvMale),JuvMale:=0]
  meso[is.na(AdMale),AdMale:=0]
  meso[is.na(EstimatedEmbryosTot),EstimatedEmbryosTot:=0]
  meso$totmale <- meso$JuvMale + meso$AdMale
  meso$EmbDividebyTot <- meso$EstimatedEmbryosTot/meso$TotalInd
  meso$EppDividebyTot <- meso$LooseEppTotal/meso$TotalInd
  meso$SCrep <- ifelse(meso$SC=="A" & meso$Clone=="D8179", "1", ifelse(meso$SC=="C" & meso$Clone=="D8222", "1", "2"))


### Restrict to weeks 1-7 and SCs A and C (not looking at AxC)

  meso7AC <- meso[Week<8 & SC=="A" | Week<8 & SC=="C"]

### Figures
  week7ACtotalind <- ggplot(data=meso7AC, aes(x=Week, y=TotalInd, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACEstEmb <- ggplot(data=meso7AC, aes(x=Week, y=EmbDividebyTot, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACPropTotMale <- ggplot(data=meso7AC, aes(x=Week, y=PropTotMale, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


### Other possible figures/panels

  week7ACTotMale <- ggplot(data=meso7AC, aes(x=Week, y=totmale, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


  week7ACEppTot <- ggplot(data=meso7AC, aes(x=Week, y=LooseEppTotal, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACEppTotbyCumPop <- ggplot(data=meso7AC, aes(x=CumInd, y=LooseEppTotal, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACTotMalebyCumPop <- ggplot(data=meso7AC, aes(x=CumInd, y=totmale, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACEmbDividebyTot <- ggplot(data=meso7AC, aes(x=Week, y=EmbDividebyTot, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACEppDividebyTot<- ggplot(data=meso7AC, aes(x=Week, y=EppDividebyTot, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACEppDividebyTotbyCumPop<- ggplot(data=meso7AC, aes(x=CumInd, y=EppDividebyTot, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACMaleDividebyTotbyCumPop<- ggplot(data=meso7AC, aes(x=CumInd, y=PropTotMale, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
