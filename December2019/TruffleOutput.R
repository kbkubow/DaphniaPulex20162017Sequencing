#!/usr/bin/env Rscript

### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(ggbeeswarm)
  library(tidyverse)

### Read in KING output file
  #kinship <- fread("polyDorsetMAF03.kin")
  kinship <- fread("polyDorsetMAF005.kin")
  #kinship <- fread("polyDorsetMAF005LDprune.kin")

### Read in SC file
  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/Superclones201617182019withObtusaandPulicaria_kingcorr_20200402_wmedrd.txt")
  #sc <- fread("CloneInfoFilePulexandObtusa_withmedrd_20200207")
  sc <- sc[Species=="pulex" & Nonindependent==0]

### Add SC info to king output
  scsubA <- data.table(ID1=sc$clone, SCA=sc$SC, medrdA=sc$medrd)
  scsubB <- data.table(ID2=sc$clone, SCB=sc$SC, medrdB=sc$medrd)
  setkey(kinship, ID2)
  setkey(scsubB, ID2)
  tmpm <- merge(kinship, scsubB)
  setkey(tmpm, ID1)
  setkey(scsubA, ID1)
  kinshipsc <- merge(tmpm, scsubA)

### Add in year and population info
  temp <- unlist(strsplit(as.character(kinshipsc$ID1), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  kinshipsc$populationA <- matdat$V3
  kinshipsc$yearA <- matdat$V2

  temp <- unlist(strsplit(as.character(kinshipsc$ID2), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  kinshipsc$populationB <- matdat$V3
  kinshipsc$yearB <- matdat$V2

  kinshipsc$pondcompare <- paste(kinshipsc$populationA, kinshipsc$populationB, sep="_")
  kinshipsc$SCcompare <- paste(kinshipsc$SCA, kinshipsc$SCB, sep="_")
  kinshipsc$clone <- ifelse(kinshipsc$SCA!="OO" & kinshipsc$SCB!="OO" &
    kinshipsc$SCA!="G" & kinshipsc$SCB!="G" &
    kinshipsc$SCA==kinshipsc$SCB, 1, 0)

  ggplot(data=kinshipsc[Kinship > 0 & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(clone))) + geom_point()
  ggplot(data=kinshipsc[Kinship > 0.2], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()

  kinshipsc$pondcompare <- paste(kinshipsc$populationA,kinshipsc$populationB, sep="_")
  kinshipsc$pondcompareB <- ifelse(kinshipsc$pondcompare=="DBunk_D8", "D8_DBunk", ifelse(
    kinshipsc$pondcompare=="DCat_D8", "D8_DCat", ifelse(kinshipsc$pondcompare=="DCat_DBunk",
    "DBunk_DCat", kinshipsc$pondcompare)
  ))

  ggplot(data=kinshipsc[Kinship > 0 & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(clone))) +
    geom_point() + facet_wrap(~pondcompareB)

    ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(clone))) +
      geom_point() + facet_wrap(~pondcompareB)

  load("./../mclonecountswide_ACfixedSNPs_20200417.Rdata")
  mclonessubID1 <- data.table(ID1=mclonecountswide$clone, ID1_ACF1=mclonecountswide$ACF1hybrid)
  mclonessubID2 <- data.table(ID2=mclonecountswide$clone, ID2_ACF1=mclonecountswide$ACF1hybrid)
  setkey(mclonessubID2, ID2)
  setkey(kinshipsc, ID2)
  mkinshipsc <- merge(kinshipsc, mclonessubID2)
  setkey(mclonessubID1, ID1)
  setkey(mkinshipsc, ID1)
  m2kinshipsc <- merge(mkinshipsc, mclonessubID1)
  m2kinshipsc$ACF1comp <- ifelse(m2kinshipsc$ID2_ACF1=="1" & m2kinshipsc$ID1_ACF1=="1" & m2kinshipsc$clone==0, 1, 0)
  m2kinshipsc$F1vsACparent <- ifelse(m2kinshipsc$SCA=="A" & m2kinshipsc$ID2_ACF1==1 |
    m2kinshipsc$SCA=="C" & m2kinshipsc$ID2_ACF1==1 | m2kinshipsc$SCB=="A" & m2kinshipsc$ID1_ACF1==1 |
    m2kinshipsc$SCB=="C" & m2kinshipsc$ID1_ACF1==1, 1, 0)

  ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(ACF1comp))) +
    geom_point() + facet_wrap(~pondcompareB)
  ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(F1vsACparent))) +
    geom_point() + facet_wrap(~pondcompareB)

  m2kinshipsc$ACrelation <- ifelse(m2kinshipsc$ACF1comp==1, "ACF1vsACF1", ifelse(
    m2kinshipsc$F1vsACparent==1, "AorCvsACF1", ifelse(m2kinshipsc$SCcompare=="A_C" |
    m2kinshipsc$SCcompare=="C_A", "AvsC", 0)))

    m2kinshipsc$ACrelation <- ifelse(m2kinshipsc$ACF1comp==1, "ACF1vsACF1", ifelse(
      m2kinshipsc$F1vsACparent==1, "AorCvsACF1", ifelse(m2kinshipsc$SCcompare=="A_C" |
      m2kinshipsc$SCcompare=="C_A", "AvsC", ifelse(m2kinshipsc$SCcompare=="C_C" |
      m2kinshipsc$SCcompare=="A_A", "AvsAorCvsC", "neither"))))


  ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(ACrelation))) +
    geom_point() + facet_wrap(~pondcompareB)

  ggplot(data=m2kinshipsc[populationA=="D8" & populationB=="D8" & medrdA>9 & medrdB>9], aes(x=IBS0, y=Kinship, color=as.factor(ACrelation))) +
    geom_point()

  ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
    aes(x=IBS0, y=Kinship, color=as.factor(ACrelation))) +
    geom_point()

    ggplot(data=m2kinshipsc[IBS0 < 0.025 & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
      aes(x=IBS0, y=Kinship, color=as.factor(ACrelation))) +
      geom_point()

      ggplot(data=m2kinshipsc[IBS0 < 0.025 & medrdA > 5 & medrdB > 5 & populationA=="D8" & populationB=="D8"],
        aes(x=IBS0, y=Kinship, color=as.factor(ACrelation))) +
        geom_point()

  m2kinshipsc$ACrelationB <- ifelse(m2kinshipsc$SCcompare=="H_C" |
    m2kinshipsc$SCcompare=="C_H", "CvsH", ifelse(m2kinshipsc$SCcompare=="W_C" |
    m2kinshipsc$SCcompare=="C_W", "CvsW", ifelse(m2kinshipsc$SCcompare=="W_H" |
    m2kinshipsc$SCcompare=="H_W", "HvsW", m2kinshipsc$ACrelation)))

    ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
      aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
      geom_point() + facet_wrap(~ACrelationB)

      ggplot(data=m2kinshipsc[Kinship > 0.2 & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
        aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
        geom_point()


        ggplot(data=m2kinshipsc[IBS0 < 0.02 & Kinship > 0.2 & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
          aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
          geom_point()

        ggplot(data=m2kinshipsc[IBS0 < 0.03 & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
          aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
          geom_point()

  m2kinshipsc$yearcompare <- paste(m2kinshipsc$yearA, m2kinshipsc$yearB, sep="_")
  m2kinshipsc$yearcompareB <- ifelse(m2kinshipsc$yearcompare=="2017_2016", "2016_2017", ifelse(
    m2kinshipsc$yearcompare=="2018_2016", "2016_2018", ifelse(
    m2kinshipsc$yearcompare=="2019_2016", "2016_2019", ifelse(
    m2kinshipsc$yearcompare=="2018_2017", "2017_2018", ifelse(
    m2kinshipsc$yearcompare=="2019_2017", "2017_2019", ifelse(
    m2kinshipsc$yearcompare=="2019_2018", "2018_2019", m2kinshipsc$yearcompare
      ))))))

        ggplot(data=m2kinshipsc[yearA==yearB & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
          aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
          geom_point() + facet_wrap(~yearcompare)

          ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
            aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
            geom_point() + facet_wrap(~yearcompareB)

            ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9 & populationA=="DBunk" & populationB=="DBunk"],
              aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
              geom_point() + facet_wrap(~yearcompareB)

              ggplot(data=m2kinshipsc[medrdA > 5 & medrdB > 5 & populationA=="DCat" & populationB=="DCat"],
                aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
                geom_point() + facet_wrap(~yearcompareB)

                ggplot(data=m2kinshipsc[medrdA > 5 & medrdB > 5 & populationA=="D8" & populationB=="D8"],
                  aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
                  geom_point() + facet_wrap(~yearcompareB)


  PO1 <- m2kinshipsc[Kinship > 0.36 & Kinship < 0.4 & IBS0 < 0.0025]
  PO2 <- m2kinshipsc[Kinship > 0.325 & Kinship < 0.36 & IBS0 < 0.0075]
  PO3 <- m2kinshipsc[Kinship > 0.24 & Kinship < 0.325 & IBS0 < 0.01]
  LikelyClones <- kinshipsc[Kinship > 0.375 & IBS0 < 0.00025]
  LikelyClonesB <- kinshipsc[Kinship > 0.35 & IBS0 < 0.00025]
  LikelyClonesMissed <- LikelyClonesB[clone!=1 & medrdA > 4 & medrdB > 4]

  PO <- rbind(PO1, PO2, PO3)


  PO1 <- kinshipsc[Kinship > 0.3 & Kinship < 0.375 & IBS0 < 0.0019]
  PO2 <- kinshipsc[Kinship > 0.25 & Kinship < 0.3 & IBS0 < 0.0020]
  PO3 <- kinshipsc[Kinship > 0.2 & Kinship < 0.25 & IBS0 < 0.0032]

  kinshipscwinpond <- kinshipsc[populationA==populationB]

  ggplot(data=kinshipscwinpond[Kinship > 0 & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(clone))) +
    geom_point() + facet_wrap(~populationA)

  kinshipscwinpond$AvsC <- ifelse(kinshipscwinpond$SCA=="A" & kinshipscwinpond$SCB=="C" |
    kinshipscwinpond$SCA=="C" & kinshipscwinpond$SCB=="A", 1, ifelse(kinshipscwinpond$ID2=="Spring_2016_D8_8.23"
    & kinshipscwinpond$ID1=="April_2017_D8_214", 2, 0))

  ggplot(data=kinshipscwinpond[Kinship > 0 & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(AvsC))) +
    geom_point() + facet_wrap(~populationA)

  ggplot(data=kinshipscwinpond[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(AvsC))) +
    geom_point() + facet_wrap(~populationA)

    ggplot(data=kinshipscwinpond[medrdA > 9 & medrdB > 9 & populationA=="D8"], aes(x=IBS0, y=Kinship, color=as.factor(AvsC))) +
      geom_point() + facet_wrap(~AvsC)

  kinshipsc$betweenF1 <- ifelse(kinshipsc$SCcompare=="O_Q" | kinshipsc$SCcompare=="Q_O", 1, ifelse(
    kinshipsc$SCcompare=="O_S" | kinshipsc$SCcompare=="S_O", 2, ifelse(
    kinshipsc$SCcompare=="Q_S" | kinshipsc$SCcompare=="S_Q", 3, ifelse(
    kinshipsc$SCcompare=="O_T" | kinshipsc$SCcompare=="T_O", 4, ifelse(
    kinshipsc$SCcompare=="S_T" | kinshipsc$SCcompare=="T_S", 5, ifelse(
    kinshipsc$SCcompare=="Q_T" | kinshipsc$SCcompare=="T_Q", 6, ifelse(
    kinshipsc$SCcompare=="O_U" | kinshipsc$SCcompare=="U_O", 7, ifelse(
    kinshipsc$SCcompare=="S_U" | kinshipsc$SCcompare=="U_S", 8, ifelse(
    kinshipsc$SCcompare=="Q_U" | kinshipsc$SCcompare=="U_Q", 9, ifelse(
    kinshipsc$SCcompare=="T_U" | kinshipsc$SCcompare=="U_T", 10, ifelse(
    kinshipsc$SCcompare=="O_X" | kinshipsc$SCcompare=="X_O", 11, ifelse(
    kinshipsc$SCcompare=="S_X" | kinshipsc$SCcompare=="X_S", 12, ifelse(
    kinshipsc$SCcompare=="Q_X" | kinshipsc$SCcompare=="X_Q", 13, ifelse(
    kinshipsc$SCcompare=="T_X" | kinshipsc$SCcompare=="X_T", 14, ifelse(
    kinshipsc$SCcompare=="U_X" | kinshipsc$SCcompare=="X_U", 15, 0
  )))))))))))))))

  ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9 & Kinship>0], aes(x=IBS0, y=Kinship, color=as.factor(betweenF1))) +
    geom_point() + facet_wrap(~betweenF1)

    ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9 & Kinship>0], aes(x=IBS0, y=Kinship, color=as.factor(betweenF1))) +
      geom_point()

    kinshipsc$betweenF1B <- ifelse(kinshipsc$betweenF1 > 0, 1, 0)
    ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9 & Kinship>0], aes(x=IBS0, y=Kinship, color=as.factor(betweenF1B))) +
      geom_point()

    ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9 & Kinship>0 & populationA=="D8" & populationB=="D8"],
     aes(x=IBS0, y=Kinship, color=as.factor(betweenF1B))) +
      geom_point()


  kinshipsc <- kinshipsc[populationA!="DOil" & populationA!="Dramp" & populationB!="DOil" & populationB!="Dramp"]
  kinshipsc$withinpond <- ifelse(kinshipsc$populationA==kinshipsc$populationB, 1, 0)
  ggplot(data=kinshipsc, aes(x=as.factor(pondcompare), y=Kinship, color=as.factor(withinpond))) + geom_beeswarm()

### From looking at a graph of Kinship by IBS0, we can see clear groupings falling out. One of these is the line
### along the bottom that appears to respond to clones (individuals from the same clonal lineage.
  # Spring_2016_D8_8.24 appears to belong to superclone H
  # Spring_2016_D8_8.23 appears to belong to superclone C (original B)
  # April_2017_DCat_10 appears to belong to superclone B

  kinshipsc$clone <- ifelse(kinshipsc$SCA=="C" & kinshipsc$ID1=="Spring_2016_D8_8.23" |
    kinshipsc$SCB=="C" & kinshipsc$ID1=="Spring_2016_D8_8.23" |
    kinshipsc$SCA=="C" & kinshipsc$ID2=="Spring_2016_D8_8.23" |
    kinshipsc$SCB=="C" & kinshipsc$ID2=="Spring_2016_D8_8.23", 2, kinshipsc$clone)

    kinshipsc$clone <- ifelse(kinshipsc$SCA=="B" & kinshipsc$ID1=="April_2017_DCat_10" |
      kinshipsc$SCB=="B" & kinshipsc$ID1=="April_2017_DCat_10" |
      kinshipsc$SCA=="B" & kinshipsc$ID2=="April_2017_DCat_10" |
      kinshipsc$SCB=="B" & kinshipsc$ID2=="April_2017_DCat_10", 3, kinshipsc$clone)

    kinshipsc$clone <- ifelse(kinshipsc$SCA=="H" & kinshipsc$ID1=="Spring_2016_D8_8.24" |
      kinshipsc$SCB=="H" & kinshipsc$ID1=="Spring_2016_D8_8.24" |
      kinshipsc$SCA=="H" & kinshipsc$ID2=="Spring_2016_D8_8.24" |
      kinshipsc$SCB=="H" & kinshipsc$ID2=="Spring_2016_D8_8.24", 4, kinshipsc$clone)



  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.19], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()
  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.3], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()

  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.3], aes(x=Kinship, y=IBS0, color=as.factor(clone))) +
    geom_point() + facet_wrap(~clone)

  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.3], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()
  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.19], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()
  ggplot(data=kinshipsc[N_SNP > 490000 & Kinship > 0.19], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()



  msc$pondcompare <- str_replace(msc$pondcompare, "D8_D10", "D10_D8")
  msc$pondcompare <- str_replace(msc$pondcompare, "DBunk_D10", "D10_DBunk")
  msc$pondcompare <- str_replace(msc$pondcompare, "DCat_D10", "D10_DCat")
  msc$pondcompare <- str_replace(msc$pondcompare, "DBunk_D8", "D8_DBunk")
  msc$pondcompare <- str_replace(msc$pondcompare, "DCat_D8", "D8_DCat")
  msc$pondcompare <- str_replace(msc$pondcompare, "W1_D8", "D8_W1")
  msc$pondcompare <- str_replace(msc$pondcompare, "W6_D8", "D8_W6")
  msc$pondcompare <- str_replace(msc$pondcompare, "W1_DBunk", "DBunk_W1")
  msc$pondcompare <- str_replace(msc$pondcompare, "W6_DBunk", "DBunk_W6")
  msc.ag <- msc[,list(meanIBS=mean(Kinship)),
    list(pondcompare, window, chr, start, stop)]
