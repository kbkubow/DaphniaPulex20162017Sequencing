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

### Let's first look at total population size over time
  ### Cleaning up and adding new variables
  meso$Clone <- as.factor(meso$Clone)
  meso$TotalPop <- meso$TotalInd*15
  meso[is.na(JuvMale),JuvMale:=0]
  meso[is.na(AdMale),AdMale:=0]
  meso[is.na(EstimatedEmbryosTot),EstimatedEmbryosTot:=0]
  meso$totmale <- meso$JuvMale + meso$AdMale
  meso$EmbDividebyTot <- meso$EstimatedEmbryosTot/meso$TotalPop
  meso$EppDividebyTot <- meso$LooseEppTotal/meso$TotalPop
  meso$SCrep <- ifelse(meso$SC=="A" & meso$Clone=="D8179", "1", ifelse(meso$SC=="C" & meso$Clone=="D8222", "1", "2"))
  meso[is.na(MomPE_5),MomPE_5:=0]
  meso[is.na(MomPE_10),MomPE_10:=0]
  meso[is.na(MomPE_15),MomPE_15:=0]
  meso[is.na(MomPE_20),MomPE_20:=0]
  meso[is.na(MomwEpp),MomwEpp:=0]
  meso[is.na(Barren),Barren:=0]
  meso$totmomPE <- meso$MomPE_5 + meso$MomPE_10 + meso$MomPE_15 + meso$MomPE_20
  meso$totalADfemale <- meso$totmomPE + meso$MomwEpp + meso$Barren
  meso$propPEtotalind <- meso$totmomPE/meso$TotalInd
  meso$propPEmoms <- meso$totmomPE/meso$totalADfemale
  meso$propMomEpp <- meso$MomwEpp/meso$TotalInd
  meso$propfill <- meso$EstimatedEmbryosTot/(meso$LooseEppTotal*2)



  meso7 <- meso[Week<8]
  meso7AC <- meso7[SC=="A" | SC=="C"]
  mesoAC <- meso[SC=="A" | SC=="C"]
  meso7AC <- meso[Week<8 & SC=="A" | Week<8 & SC=="C"]

  week7ACtotalind <- ggplot(data=meso7AC, aes(x=Week, y=TotalPop, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("# Individuals")

  t1 <- lmer(TotalPop~1+(1|Week)+ (1|Clone)+SC, data=meso7AC)
  t2 <- lmer(TotalPop~1+(1|Week)+ (1|Clone), data=meso7AC)

  anova(t1, t2)

  t1 <- lmer(TotalPop~1+(1|Clone)+SC, data=meso7AC[Week==7])
  t2 <- lmer(TotalPop~1+(1|Clone), data=meso7AC[Week==7])

  anova(t1, t2)

  week7ACmomPEnormWeek <- ggplot(data=meso7AC, aes(x=Week, y=propPEtotalind, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Asexual Females")

  t1 <- glmer(propPEtotalind~1+(1|Week)+ (1|Clone)+SC, data=meso7AC, weights=TotalInd, family=binomial)
  t2 <- glmer(propPEtotalind~1+(1|Week)+ (1|Clone), data=meso7AC, weights=TotalInd, family=binomial)

  anova(t1, t2)

  week7ACMomwEppnormWeek <- ggplot(data=meso7AC, aes(x=Week, y=propMomEpp, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Sexual Females")

  t1 <- glmer(propMomEpp~1+(1|Week)+ (1|Clone)+SC, data=meso7AC, weights=TotalInd, family=binomial)
  t2 <- glmer(propMomEpp~1+(1|Week)+ (1|Clone), data=meso7AC, weights=TotalInd, family=binomial)

  anova(t1, t2)

  week7ACPropMalenormWeek <- ggplot(data=meso7AC, aes(x=Week, y=PropTotMale, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Male")

  t1 <- glmer(PropTotMale~1+(1|Week)+ (1|Clone)+SC, data=meso7AC, weights=TotalInd, family=binomial)
  t2 <- glmer(PropTotMale~1+(1|Week)+ (1|Clone), data=meso7AC, weights=TotalInd, family=binomial)

  anova(t1, t2)

  week7ACEmbWeek <- ggplot(data=meso7AC, aes(x=Week, y=log10(EstimatedEmbryosTot), color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("# Sexual Embryos")

  t1 <- lmer(EstimatedEmbryosTot~1+(1|Week)+ (1|Clone) + SC, data=meso7AC)
  t2 <- lmer(EstimatedEmbryosTot~1+(1|Week)+ (1|Clone), data=meso7AC)

  anova(t1, t2)
