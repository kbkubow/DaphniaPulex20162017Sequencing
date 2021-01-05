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

            load("mscACvsrest.ag_G.Rdata")
            totalibsag <- mscACvsrest.ag

            load("mscACvsrest_POBC_G.Rdata")
            totalibsagPO <- mscACvsrest_POBC

            totalibsagPO$ACvspond <- ifelse(totalibsagPO$ACvspond=="A_obtusa" & totalibsagPO$SCA_C=="C", "C_obtusa", ifelse(
              totalibsagPO$ACvspond=="A_pulicaria" & totalibsagPO$SCA_C=="C", "C_pulicaria", totalibsagPO$ACvspond
            ))

            load("D8DBunkDCat.ag_G.Rdata")
            totalibsagD8DBunk <- D8DBunkDCat.ag

            load("D8DBunkDCatvsD10W.ag_G.Rdata")
            totalibsD10W <- D8DBunkDCatvsD10W.ag

            load("D8DBunkDCatvsPO.ag_G.Rdata")
            totalibsObPul <- D8DBunkDCatvsPO.ag

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

  AvsC <- total[ACvspondB=="A_D8DBunk_C_D8DBunk"]
  AvsCsub <- AvsC[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
  colnames(AvsCsub) <- c("window", "chr", "start", "stop", "AvsCmeanIBS")

  D8_D8 <- total[ACvspondB=="D8_D8"]
  D8_D8$SCcompare <- paste(D8_D8$SCA_C, D8_D8$SCB_C, sep="_vs_")
  D8D8SCcomp <- unique(D8_D8$SCcompare)

  AvsC_D8D8 <- foreach(i=1:length(D8D8SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8D8SCcomp[i]
    tmp <- D8_D8[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_D8D8", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  DBunk_DBunk <- total[ACvspondB=="DBunk_DBunk"]
  DBunk_DBunk$SCcompare <- paste(DBunk_DBunk$SCA_C, DBunk_DBunk$SCB_C, sep="_vs_")
  DBunkDBunkSCcomp <- unique(DBunk_DBunk$SCcompare)

  AvsC_DBunkDBunk <- foreach(i=1:length(DBunkDBunkSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- DBunkDBunkSCcomp[i]
    tmp <- DBunk_DBunk[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_DBunkDBunk", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  DCat_DCat <- total[ACvspondB=="DCat_DCat"]
  DCat_DCat$SCcompare <- paste(DCat_DCat$SCA_C, DCat_DCat$SCB_C, sep="_vs_")
  DCatDCatSCcomp <- unique(DCat_DCat$SCcompare)

  AvsC_DCatDCat <- foreach(i=1:length(DCatDCatSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- DCatDCatSCcomp[i]
    tmp <- DCat_DCat[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_DCatDCat", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  D8_DBunk <- total[ACvspondB=="D8_DBunk"]
  D8_DBunk$SCcompare <- paste(D8_DBunk$SCA_C, D8_DBunk$SCB_C, sep="_vs_")
  D8DBunkSCcomp <- unique(D8_DBunk$SCcompare)

  AvsC_D8DBunk <- foreach(i=1:length(D8DBunkSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DBunkSCcomp[i]
    tmp <- D8_DBunk[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_D8DBunk", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  D8_DCat <- total[ACvspondB=="D8_DCat"]
  D8_DCat$SCcompare <- paste(D8_DCat$SCA_C, D8_DCat$SCB_C, sep="_vs_")
  D8DCatSCcomp <- unique(D8_DCat$SCcompare)

  AvsC_D8DCat <- foreach(i=1:length(D8DCatSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DCatSCcomp[i]
    tmp <- D8_DCat[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_D8DCat", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  DBunk_DCat <- total[ACvspondB=="DBunk_DCat"]
  DBunk_DCat$SCcompare <- paste(DBunk_DCat$SCA_C, DBunk_DCat$SCB_C, sep="_vs_")
  DBunkDCatSCcomp <- unique(DBunk_DCat$SCcompare)

  AvsC_DBunkDCat <- foreach(i=1:length(DBunkDCatSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- DBunkDCatSCcomp[i]
    tmp <- DBunk_DCat[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_DBunkDCat", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  totalttestACotherpaired <- rbind(AvsC_D8D8, AvsC_DBunkDBunk, AvsC_DCatDCat, AvsC_D8DBunk, AvsC_D8DCat, AvsC_DBunkDCat)

  save(totalttestACotherpaired, file="totalttestACotherpaired_20200817.Rdata")


  A_D10 <- total[ACvspondB=="A_D10"]
  A_D10$SCcompare <- paste(A_D10$SCA_C, A_D10$SCB_C, sep="_vs_")
  A_D10SCcomp <- unique(A_D10$SCcompare)

  AvsC_A_D10 <- foreach(i=1:length(A_D10SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- A_D10SCcomp[i]
    tmp <- A_D10[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_A_D10", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  A_W1 <- total[ACvspondB=="A_W1"]
  A_W1$SCcompare <- paste(A_W1$SCA_C, A_W1$SCB_C, sep="_vs_")
  A_W1SCcomp <- unique(A_W1$SCcompare)

  AvsC_A_W1 <- foreach(i=1:length(A_W1SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- A_W1SCcomp[i]
    tmp <- A_W1[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_A_W1", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  A_W6 <- total[ACvspondB=="A_W6"]
  A_W6$SCcompare <- paste(A_W6$SCA_C, A_W6$SCB_C, sep="_vs_")
  A_W6SCcomp <- unique(A_W6$SCcompare)

  AvsC_A_W6 <- foreach(i=1:length(A_W6SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- A_W6SCcomp[i]
    tmp <- A_W6[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_A_W6", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  A_pulicaria <- total[ACvspondB=="A_pulicaria"]
  A_pulicaria$SCcompare <- paste(A_pulicaria$SCA_C, A_pulicaria$SCB_C, sep="_vs_")
  A_pulicariaSCcomp <- unique(A_pulicaria$SCcompare)

  AvsC_A_pulicaria <- foreach(i=1:length(A_pulicariaSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- A_pulicariaSCcomp[i]
    tmp <- A_pulicaria[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_A_pulicaria", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  A_obtusa <- total[ACvspondB=="A_obtusa"]
  A_obtusa$SCcompare <- paste(A_obtusa$SCA_C, A_obtusa$SCB_C, sep="_vs_")
  A_obtusaSCcomp <- unique(A_obtusa$SCcompare)

  AvsC_A_obtusa <- foreach(i=1:length(A_obtusaSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- A_obtusaSCcomp[i]
    tmp <- A_obtusa[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_A_obtusa", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  C_D10 <- total[ACvspondB=="C_D10"]
  C_D10$SCcompare <- paste(C_D10$SCA_C, C_D10$SCB_C, sep="_vs_")
  C_D10SCcomp <- unique(C_D10$SCcompare)

  AvsC_C_D10 <- foreach(i=1:length(C_D10SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- C_D10SCcomp[i]
    tmp <- C_D10[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_C_D10", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  C_W1 <- total[ACvspondB=="C_W1"]
  C_W1$SCcompare <- paste(C_W1$SCA_C, C_W1$SCB_C, sep="_vs_")
  C_W1SCcomp <- unique(C_W1$SCcompare)

  AvsC_C_W1 <- foreach(i=1:length(C_W1SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- C_W1SCcomp[i]
    tmp <- C_W1[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_C_W1", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  C_W6 <- total[ACvspondB=="C_W6"]
  C_W6$SCcompare <- paste(C_W6$SCA_C, C_W6$SCB_C, sep="_vs_")
  C_W6SCcomp <- unique(C_W6$SCcompare)

  AvsC_C_W6 <- foreach(i=1:length(C_W6SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- C_W6SCcomp[i]
    tmp <- C_W6[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_C_W6", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  C_pulicaria <- total[ACvspondB=="C_pulicaria"]
  C_pulicaria$SCcompare <- paste(C_pulicaria$SCA_C, C_pulicaria$SCB_C, sep="_vs_")
  C_pulicariaSCcomp <- unique(C_pulicaria$SCcompare)

  AvsC_C_pulicaria <- foreach(i=1:length(C_pulicariaSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- C_pulicariaSCcomp[i]
    tmp <- C_pulicaria[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_C_pulicaria", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  C_obtusa <- total[ACvspondB=="C_obtusa"]
  C_obtusa$SCcompare <- paste(C_obtusa$SCA_C, C_obtusa$SCB_C, sep="_vs_")
  C_obtusaSCcomp <- unique(C_obtusa$SCcompare)

  AvsC_C_obtusa <- foreach(i=1:length(C_obtusaSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- C_obtusaSCcomp[i]
    tmp <- C_obtusa[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
    stats <- data.table(comppond="AvsC_C_obtusa", compsc=sc, t=ttest$statistic, df=ttest$parameter,
    p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
    meandiff=ttest$estimate[1])
  }

  totalttestACACoutgrouppaired <- rbind(AvsC_A_D10, AvsC_A_W1, AvsC_A_W6, AvsC_A_pulicaria, AvsC_A_obtusa,
    AvsC_C_D10, AvsC_C_W1, AvsC_C_W6, AvsC_C_pulicaria, AvsC_C_obtusa)

  save(totalttestACACoutgrouppaired, file="totalttestACACoutgrouppaired_20200817.Rdata")


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

              load("mscACvsrest.ag_G.Rdata")
              totalibsag <- mscACvsrest.ag

              load("mscACvsrest_POBC_G.Rdata")
              totalibsagPO <- mscACvsrest_POBC

              totalibsagPO$ACvspond <- ifelse(totalibsagPO$ACvspond=="A_obtusa" & totalibsagPO$SCA_C=="C", "C_obtusa", ifelse(
                totalibsagPO$ACvspond=="A_pulicaria" & totalibsagPO$SCA_C=="C", "C_pulicaria", totalibsagPO$ACvspond
              ))

              load("D8DBunkDCat.ag_G.Rdata")
              totalibsagD8DBunk <- D8DBunkDCat.ag

              load("D8DBunkDCatvsD10W.ag_G.Rdata")
              totalibsD10W <- D8DBunkDCatvsD10W.ag

              load("D8DBunkDCatvsPO.ag_G.Rdata")
              totalibsObPul <- D8DBunkDCatvsPO.ag

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

    AvsC <- total[ACvspondB=="A_D8DBunk_C_D8DBunk"]
    AvsCsub <- AvsC[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(AvsCsub) <- c("window", "chr", "start", "stop", "AvsCmeanIBS")

    D8DBunkDCatvsD10 <- total[ACvspondB=="D8DBunkDCatvsD10"]
    D8DBunkDCatvsD10$SCcompare <- paste(D8DBunkDCatvsD10$SCA_C, D8DBunkDCatvsD10$SCB_C, sep="_vs_")
    D8DBunkDCatvsD10SCcomp <- unique(D8DBunkDCatvsD10$SCcompare)

    AvsC_D8DBunkDCatvsD10 <- foreach(i=1:length(D8DBunkDCatvsD10SCcomp), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      sc <- D8DBunkDCatvsD10SCcomp[i]
      tmp <- D8DBunkDCatvsD10[SCcompare==sc]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsD10", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meandiff=ttest$estimate[1])
    }

    D8DBunkDCatvsW1 <- total[ACvspondB=="D8DBunkDCatvsW1"]
    D8DBunkDCatvsW1$SCcompare <- paste(D8DBunkDCatvsW1$SCA_C, D8DBunkDCatvsW1$SCB_C, sep="_vs_")
    D8DBunkDCatvsW1SCcomp <- unique(D8DBunkDCatvsW1$SCcompare)

    AvsC_D8DBunkDCatvsW1 <- foreach(i=1:length(D8DBunkDCatvsW1SCcomp), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      sc <- D8DBunkDCatvsW1SCcomp[i]
      tmp <- D8DBunkDCatvsW1[SCcompare==sc]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsW1", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meandiff=ttest$estimate[1])
    }

    D8DBunkDCatvsW6 <- total[ACvspondB=="D8DBunkDCatvsW6"]
    D8DBunkDCatvsW6$SCcompare <- paste(D8DBunkDCatvsW6$SCA_C, D8DBunkDCatvsW6$SCB_C, sep="_vs_")
    D8DBunkDCatvsW6SCcomp <- unique(D8DBunkDCatvsW6$SCcompare)

    AvsC_D8DBunkDCatvsW1 <- foreach(i=1:length(D8DBunkDCatvsW6SCcomp), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      sc <- D8DBunkDCatvsW6SCcomp[i]
      tmp <- D8DBunkDCatvsW6[SCcompare==sc]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsW6", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meandiff=ttest$estimate[1])
    }

    D8DBunkDCatvsW6 <- total[ACvspondB=="D8DBunkDCatvsW6"]
    D8DBunkDCatvsW6$SCcompare <- paste(D8DBunkDCatvsW6$SCA_C, D8DBunkDCatvsW6$SCB_C, sep="_vs_")
    D8DBunkDCatvsW6SCcomp <- unique(D8DBunkDCatvsW6$SCcompare)

    AvsC_D8DBunkDCatvsW6 <- foreach(i=1:length(D8DBunkDCatvsW6SCcomp), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      sc <- D8DBunkDCatvsW6SCcomp[i]
      tmp <- D8DBunkDCatvsW6[SCcompare==sc]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsW6", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meandiff=ttest$estimate[1])
    }

    D8DBunkDCatvsPulicaria <- total[ACvspondB=="D8DBunkDCatvsPulicaria"]
    D8DBunkDCatvsPulicaria$SCcompare <- paste(D8DBunkDCatvsPulicaria$SCA_C, D8DBunkDCatvsPulicaria$SCB_C, sep="_vs_")
    D8DBunkDCatvsPulicariaSCcomp <- unique(D8DBunkDCatvsPulicaria$SCcompare)

    AvsC_D8DBunkDCatvsPulicaria <- foreach(i=1:length(D8DBunkDCatvsPulicariaSCcomp), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      sc <- D8DBunkDCatvsPulicariaSCcomp[i]
      tmp <- D8DBunkDCatvsPulicaria[SCcompare==sc]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      ttest <- t.test(m$AvsCmeanIBS, m$tmp, paired=TRUE)
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsPulicaria", compsc=sc, t=ttest$statistic, df=ttest$parameter,
      p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
      meandiff=ttest$estimate[1])
    }


  load("totalttest_20200812.Rdata")

  ggplot(data=totalttest, aes(x=t, fill=as.factor(sig0.5))) + geom_histogram(binwidth=5) +
  facet_wrap(~comppond, scales="free_y") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
