#!/usr/bin/env Rscript

    ### libraries
      library(data.table)
      library(foreach)
      library(ggplot2)
      library(tidyverse)
      library(doMC)
      registerDoMC(20)


  ### Load processed IBS by window files
      inputfiles <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="mscACvsrest.ag")

      totalibsag <- foreach(i=1:length(inputfiles), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfiles[i], sep="")
        load(f)
        mscACvsrest.ag
      }

      inputfilesPO <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="mscACvsrest_POBC")

      totalibsagPO <- foreach(i=1:length(inputfilesPO), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesPO[i], sep="")
        load(f)
        mscACvsrest_POBC
      }

      totalibsagPO$ACvspond <- ifelse(totalibsagPO$ACvspond=="A_obtusa" & totalibsagPO$SCA_C=="C", "C_obtusa", ifelse(
        totalibsagPO$ACvspond=="A_pulicaria" & totalibsagPO$SCA_C=="C", "C_pulicaria", totalibsagPO$ACvspond
      ))

      inputfilesDBunkD8 <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="D8DBunkDCat.ag")

      totalibsagD8DBunk <- foreach(i=1:length(inputfilesDBunkD8), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesDBunkD8[i], sep="")
        load(f)
        D8DBunkDCat.ag
      }


      inputfilesObtusaPulicaria <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="D8DBunkDCatvsPO.ag")

      totalibsObPul <- foreach(i=1:length(inputfilesObtusaPulicaria), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesObtusaPulicaria[i], sep="")
        load(f)
        D8DBunkDCatvsPO.ag
      }

      inputfilesD10W <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="D8DBunkDCatvsD10W.ag")

      totalibsD10W <- foreach(i=1:length(inputfilesD10W), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesD10W[i], sep="")
        load(f)
        D8DBunkDCatvsD10W.ag
      }

      ### Read in only G files first to figure out how to make the figure

            #load("mscACvsrest.ag_G.Rdata")
            #totalibsag <- mscACvsrest.ag

            #load("mscACvsrest_POBC_G.Rdata")
            #totalibsagPO <- mscACvsrest_POBC

            #totalibsagPO$ACvspond <- ifelse(totalibsagPO$ACvspond=="A_obtusa" & totalibsagPO$SCA_C=="C", "C_obtusa", ifelse(
            #  totalibsagPO$ACvspond=="A_pulicaria" & totalibsagPO$SCA_C=="C", "C_pulicaria", totalibsagPO$ACvspond
            #))

            #load("D8DBunkDCat.ag_G.Rdata")
            #totalibsagD8DBunk <- D8DBunkDCat.ag

            #load("D8DBunkDCatvsD10W.ag_G.Rdata")
            #totalibsD10W <- D8DBunkDCatvsD10W.ag

            #load("D8DBunkDCatvsPO.ag_G.Rdata")
            #totalibsObPul <- D8DBunkDCatvsPO.ag

            total <- rbind(totalibsag, totalibsagPO, totalibsagD8DBunk, totalibsD10W, totalibsObPul)

            total$ACvspondB <- ifelse(total$ACvspond=="D8_D10" | total$ACvspond=="DBunk_D10" |
              total$ACvspond=="DCat_D10", "D8DBunkDCatvsD10", ifelse(total$ACvspond=="D8_W1" |
              total$ACvspond=="DBunk_W1" | total$ACvspond=="DCat_W1", "D8DBunkDCatvsW1", ifelse(
              total$ACvspond=="D8_W6" | total$ACvspond=="DBunk_W6" | total$ACvspond=="DCat_W6",
              "D8DBunkDCatvsW6", total$ACvspond
              )))

              total$type <- ifelse(total$ACvspondB=="A_D8DBunk_C_D8DBunk", "AvsC", ifelse(
                total$ACvspondB=="D8_D8" | total$ACvspondB=="DBunk_DBunk" | total$ACvspondB=="DCat_DCat",
                "withinpond", ifelse(total$ACvspondB=="D8_DBunk" | total$ACvspondB=="D8_DCat" |
                total$ACvspondB=="DBunk_DCat", "betweenpond", ifelse(total$ACvspondB=="A_D10" |
                total$ACvspondB=="A_W1" | total$ACvspondB=="A_W6" | total$ACvspondB=="A_pulicaria" |
                total$ACvspondB=="A_obtusa", "Avsout", ifelse(total$ACvspondB=="C_D10" |
                total$ACvspondB=="C_W1" | total$ACvspondB=="C_W6" | total$ACvspondB=="C_pulicaria" |
                total$ACvspondB=="C_obtusa", "Cvsout", "othervsout"
              )))))

              total$type <- factor(total$type, levels=c("AvsC", "withinpond", "betweenpond",
                "Avsout", "Cvsout", "othervsout"))
              total$ACvspondB <- factor(total$ACvspondB, levels=c("A_D8DBunk_C_D8DBunk",
                "D8_D8", "DBunk_DBunk", "DCat_DCat", "D8_DBunk", "D8_DCat", "DBunk_DCat", "A_D10",
                "A_W1", "A_W6", "A_pulicaria", "A_obtusa", "C_D10", "C_W1", "C_W6", "C_pulicaria",
                "C_obtusa", "D8DBunkDCatvsD10", "D8DBunkDCatvsW1", "D8DBunkDCatvsW6",
                "D8DBunkDCatvsPulicaria", "D8DBunkDCatvsObtusa"))

  AvsC <- total$meanIBS[total$ACvspondB=="A_D8DBunk_C_D8DBunk"]

  D8_D8 <- total[ACvspondB=="D8_D8"]
  D8_D8$SCcompare <- paste(D8_D8$SCA_C, D8_D8$SCB_C, sep="_vs_")
  D8D8SCcomp <- unique(D8_D8$SCcompare)

  AvsC_D8D8 <- foreach(i=1:length(D8D8SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8D8SCcomp[i]
    y <- D8_D8$meanIBS[D8_D8$SCcompare==sc]
    ttest <- t.test(AvsC, y)
    stats <- data.table(comppond="AvsC_D8D8", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meanx=ttest$estimate[1], meany=ttest$estimate[2])
  }

  DBunk_DBunk <- total[ACvspondB=="DBunk_DBunk"]
  DBunk_DBunk$SCcompare <- paste(DBunk_DBunk$SCA_C, DBunk_DBunk$SCB_C, sep="_vs_")
  DBunkDBunkSCcomp <- unique(DBunk_DBunk$SCcompare)

  AvsC_DBunkDBunk <- foreach(i=1:length(DBunkDBunkSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- DBunkDBunkSCcomp[i]
    y <- DBunk_DBunk$meanIBS[DBunk_DBunk$SCcompare==sc]
    ttest <- t.test(AvsC, y)
    stats <- data.table(comppond="AvsC_DBunkDBunk", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meanx=ttest$estimate[1], meany=ttest$estimate[2])
  }

  DCat_DCat <- total[ACvspondB=="DCat_DCat"]
  DCat_DCat$SCcompare <- paste(DCat_DCat$SCA_C, DCat_DCat$SCB_C, sep="_vs_")
  DCatDCatSCcomp <- unique(DCat_DCat$SCcompare)

  AvsC_DCatDCat <- foreach(i=1:length(DCatDCatSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- DCatDCatSCcomp[i]
    y <- DCat_DCat$meanIBS[DCat_DCat$SCcompare==sc]
    ttest <- t.test(AvsC, y)
    stats <- data.table(comppond="AvsC_DCatDCat", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meanx=ttest$estimate[1], meany=ttest$estimate[2])
  }

  D8_DBunk <- total[ACvspondB=="D8_DBunk"]
  D8_DBunk$SCcompare <- paste(D8_DBunk$SCA_C, D8_DBunk$SCB_C, sep="_vs_")
  D8DBunkSCcomp <- unique(D8_DBunk$SCcompare)

  AvsC_D8DBunk <- foreach(i=1:length(D8DBunkSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DBunkSCcomp[i]
    y <- D8_DBunk$meanIBS[D8_DBunk$SCcompare==sc]
    ttest <- t.test(AvsC, y)
    stats <- data.table(comppond="AvsC_D8DBunk", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meanx=ttest$estimate[1], meany=ttest$estimate[2])
  }

  D8_DCat <- total[ACvspondB=="D8_DCat"]
  D8_DCat$SCcompare <- paste(D8_DCat$SCA_C, D8_DCat$SCB_C, sep="_vs_")
  D8DCatSCcomp <- unique(D8_DCat$SCcompare)

  AvsC_D8DCat <- foreach(i=1:length(D8DCatSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DCatSCcomp[i]
    y <- D8_DCat$meanIBS[D8_DCat$SCcompare==sc]
    ttest <- t.test(AvsC, y)
    stats <- data.table(comppond="AvsC_D8DCat", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meanx=ttest$estimate[1], meany=ttest$estimate[2])
  }

  DBunk_DCat <- total[ACvspondB=="DBunk_DCat"]
  DBunk_DCat$SCcompare <- paste(DBunk_DCat$SCA_C, DBunk_DCat$SCB_C, sep="_vs_")
  DBunkDCatSCcomp <- unique(DBunk_DCat$SCcompare)

  AvsC_DBunkDCat <- foreach(i=1:length(DBunkDCatSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- DBunkDCatSCcomp[i]
    y <- DBunk_DCat$meanIBS[DBunk_DCat$SCcompare==sc]
    ttest <- t.test(AvsC, y)
    stats <- data.table(comppond="AvsC_DBunkDCat", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meanx=ttest$estimate[1], meany=ttest$estimate[2])
  }

  totalttest <- rbind(AvsC_D8D8, AvsC_DBunkDBunk, AvsC_DCatDCat, AvsC_D8DBunk, AvsC_D8DCat, AvsC_DBunkDCat)

  save(totalttest, file="totalttest_20200812.Rdata")

  load("totalttest_20200812.Rdata")

  ggplot(data=totalttest, aes(x=t, fill=as.factor(sig0.5))) + geom_histogram(binwidth=5) +
  facet_wrap(~comppond, scales="free_y") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
