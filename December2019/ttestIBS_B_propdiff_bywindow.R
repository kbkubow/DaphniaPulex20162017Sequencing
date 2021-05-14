#!/usr/bin/env Rscript

    ### libraries
      library(data.table)
      library(foreach)
      library(ggplot2)
      library(tidyverse)
      library(doMC)
      registerDoMC(20)


  ### Load processed IBS by window files
      #inputfiles <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="mscACvsrest.ag")
      inputfiles <- list.files(path="/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/slidingwindowIBS/", pattern="mscACvsrest.ag")

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
  D8D8window <- unique(D8_D8$window)

  AvsC_D8D8 <- foreach(i=1:length(D8D8window), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- D8D8window[i]
    tmp <- D8_D8[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8D8", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  DBunk_DBunk <- total[ACvspondB=="DBunk_DBunk"]
  DBunkDBunkwindow <- unique(DBunk_DBunk$window)

  AvsC_DBunkDBunk <- foreach(i=1:length(DBunkDBunkwindow), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- DBunkDBunkwindow[i]
    tmp <- DBunk_DBunk[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_DBunkDBunk", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  DCat_DCat <- total[ACvspondB=="DCat_DCat"]
  DCatDCatwindow <- unique(DCat_DCat$window)

  AvsC_DCatDCat <- foreach(i=1:length(DCatDCatwindow), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- DCatDCatwindow[i]
    tmp <- DCat_DCat[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_DCatDCat", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  D8_DBunk <- total[ACvspondB=="D8_DBunk"]
  D8DBunkwindow <- unique(D8_DBunk$window)

  AvsC_D8DBunk <- foreach(i=1:length(D8DBunkwindow), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- D8DBunkwindow[i]
    tmp <- D8_DBunk[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8DBunk", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  D8_DCat <- total[ACvspondB=="D8_DCat"]
  D8DCatwindow <- unique(D8_DCat$window)

  AvsC_D8DCat <- foreach(i=1:length(D8DCatwindow), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- D8DCatwindow[i]
    tmp <- D8_DCat[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8DCat", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  DBunk_DCat <- total[ACvspondB=="DBunk_DCat"]
  DBunkDCatwindow <- unique(DBunk_DCat$window)

  AvsC_DBunkDCat <- foreach(i=1:length(DBunkDCatwindow), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- DBunkDCatwindow[i]
    tmp <- DBunk_DCat[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_DBunkDCat", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  totalpropsimACotherbywindow <- rbind(AvsC_D8D8, AvsC_DBunkDBunk, AvsC_DCatDCat, AvsC_D8DBunk, AvsC_D8DCat, AvsC_DBunkDCat)

  save(totalpropsimACotherbywindow, file="totalpropsimACotherbywindow_20200824.Rdata")


  A_D10 <- total[ACvspondB=="A_D10"]
  A_D10window <- unique(A_D10$window)

  AvsC_A_D10 <- foreach(i=1:length(A_D10window), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- A_D10window[i]
    tmp <- A_D10[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_A_D10", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  A_W1 <- total[ACvspondB=="A_W1"]
  A_W1window <- unique(A_W1$window)

  AvsC_A_W1 <- foreach(i=1:length(A_W1window), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- A_W1window[i]
    tmp <- A_W1[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_A_W1", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  A_W6 <- total[ACvspondB=="A_W6"]
  A_W6window <- unique(A_W6$window)

  AvsC_A_W6 <- foreach(i=1:length(A_W6window), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- A_W6window[i]
    tmp <- A_W6[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_A_W6", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  A_pulicaria <- total[ACvspondB=="A_pulicaria"]
  A_pulicariawindow <- unique(A_pulicaria$window)

  AvsC_A_pulicaria <- foreach(i=1:length(A_pulicariawindow), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- A_pulicariawindow[i]
    tmp <- A_pulicaria[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_A_pulicaria", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  A_obtusa <- total[ACvspondB=="A_obtusa"]
  A_obtusawindow <- unique(A_obtusa$window)

  AvsC_A_obtusa <- foreach(i=1:length(A_obtusawindow), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- A_obtusawindow[i]
    tmp <- A_obtusa[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_A_obtusa", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  C_D10 <- total[ACvspondB=="C_D10"]
  C_D10window <- unique(C_D10$window)

  AvsC_C_D10 <- foreach(i=1:length(C_D10window), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- C_D10window[i]
    tmp <- C_D10[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_C_D10", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  C_W1 <- total[ACvspondB=="C_W1"]
  C_W1window <- unique(C_W1$window)

  AvsC_C_W1 <- foreach(i=1:length(C_W1window), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- C_W1window[i]
    tmp <- C_W1[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_C_W1", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  C_W6 <- total[ACvspondB=="C_W6"]
  C_W6window <- unique(C_W6$window)

  AvsC_C_W6 <- foreach(i=1:length(C_W6window), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- C_W6window[i]
    tmp <- C_W6[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_C_W6", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  C_pulicaria <- total[ACvspondB=="C_pulicaria"]
  C_pulicariawindow <- unique(C_pulicaria$window)

  AvsC_C_pulicaria <- foreach(i=1:length(C_pulicariawindow), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- C_pulicariawindow[i]
    tmp <- C_pulicaria[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_C_pulicaria", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  C_obtusa <- total[ACvspondB=="C_obtusa"]
  C_obtusawindow <- unique(C_obtusa$window)

  AvsC_C_obtusa <- foreach(i=1:length(C_obtusawindow), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    win <- C_obtusawindow[i]
    tmp <- C_obtusa[window==win]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_C_obtusa", window=win, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  totalpropsimACACoutgroupbywindow <- rbind(AvsC_A_D10, AvsC_A_W1, AvsC_A_W6, AvsC_A_pulicaria, AvsC_A_obtusa,
    AvsC_C_D10, AvsC_C_W1, AvsC_C_W6, AvsC_C_pulicaria, AvsC_C_obtusa)

  save(totalpropsimACACoutgroupbywindow, file="totalpropsimACACoutgroupbywindow_20200824.Rdata")


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

            #  load("mscACvsrest.ag_G.Rdata")
            #  totalibsag <- mscACvsrest.ag

            #  load("mscACvsrest_POBC_G.Rdata")
            #  totalibsagPO <- mscACvsrest_POBC

            #  totalibsagPO$ACvspond <- ifelse(totalibsagPO$ACvspond=="A_obtusa" & totalibsagPO$SCA_C=="C", "C_obtusa", ifelse(
            #    totalibsagPO$ACvspond=="A_pulicaria" & totalibsagPO$SCA_C=="C", "C_pulicaria", totalibsagPO$ACvspond
            #  ))

            #  load("D8DBunkDCat.ag_G.Rdata")
            #  totalibsagD8DBunk <- D8DBunkDCat.ag

            #  load("D8DBunkDCatvsD10W.ag_G.Rdata")
            #  totalibsD10W <- D8DBunkDCatvsD10W.ag

            #  load("D8DBunkDCatvsPO.ag_G.Rdata")
            #  totalibsObPul <- D8DBunkDCatvsPO.ag

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

    windowinfo <- AvsC[, c("window", "chr", "start", "stop"), with=TRUE]
    save(windowinfo, file="windowinfo_20200824.Rdata")


    D8DBunkDCatvsD10 <- total[ACvspondB=="D8DBunkDCatvsD10"]
    D8DBunkDCatvsD10window <- unique(D8DBunkDCatvsD10$window)

    AvsC_D8DBunkDCatvsD10 <- foreach(i=1:length(D8DBunkDCatvsD10window), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      win <- D8DBunkDCatvsD10window[i]
      tmp <- D8DBunkDCatvsD10[window==win]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      m$ratiodiff <- m$AvsCmeanIBS/m$tmp
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsD10", window=win, mean_ratiosim=mean(m$ratiodiff),
        sd_ratiosim=sd(m$ratiodiff))
    }

    D8DBunkDCatvsW1 <- total[ACvspondB=="D8DBunkDCatvsW1"]
    D8DBunkDCatvsW1window <- unique(D8DBunkDCatvsW1$window)

    AvsC_D8DBunkDCatvsW1 <- foreach(i=1:length(D8DBunkDCatvsW1window), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      win <- D8DBunkDCatvsW1window[i]
      tmp <- D8DBunkDCatvsW1[window==win]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      m$ratiodiff <- m$AvsCmeanIBS/m$tmp
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsW1", window=win, mean_ratiosim=mean(m$ratiodiff),
        sd_ratiosim=sd(m$ratiodiff))
    }

    D8DBunkDCatvsW6 <- total[ACvspondB=="D8DBunkDCatvsW6"]
    D8DBunkDCatvsW6window <- unique(D8DBunkDCatvsW6$window)

    AvsC_D8DBunkDCatvsW6 <- foreach(i=1:length(D8DBunkDCatvsW6window), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      win <- D8DBunkDCatvsW6window[i]
      tmp <- D8DBunkDCatvsW6[window==win]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      m$ratiodiff <- m$AvsCmeanIBS/m$tmp
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsW6", window=win, mean_ratiosim=mean(m$ratiodiff),
        sd_ratiosim=sd(m$ratiodiff))
    }

    D8DBunkDCatvsPulicaria <- total[ACvspondB=="D8DBunkDCatvsPulicaria"]
    D8DBunkDCatvsPulicariawindow <- unique(D8DBunkDCatvsPulicaria$window)

    AvsC_D8DBunkDCatvsPulicaria <- foreach(i=1:length(D8DBunkDCatvsPulicariawindow), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      win <- D8DBunkDCatvsPulicariawindow[i]
      tmp <- D8DBunkDCatvsPulicaria[window==win]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      m$ratiodiff <- m$AvsCmeanIBS/m$tmp
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsPulicaria", window=win, mean_ratiosim=mean(m$ratiodiff),
        sd_ratiosim=sd(m$ratiodiff))
    }

    D8DBunkDCatvsObtusa <- total[ACvspondB=="D8DBunkDCatvsObtusa"]
    D8DBunkDCatvsObtusawindow <- unique(D8DBunkDCatvsObtusa$window)

    AvsC_D8DBunkDCatvsObtusa <- foreach(i=1:length(D8DBunkDCatvsObtusawindow), .combine="rbind")%do%{
      print(paste("Getting data: ", i, sep=""))
      win <- D8DBunkDCatvsObtusawindow[i]
      tmp <- D8DBunkDCatvsObtusa[window==win]
      tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
      colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
      setkey(AvsCsub, window, chr, start, stop)
      setkey(tmpsub, window, chr, start, stop)
      m <- merge(AvsCsub, tmpsub)
      m$ratiodiff <- m$AvsCmeanIBS/m$tmp
      stats <- data.table(comppond="AvsC_D8DBunkDCatvsObtusa", window=win, mean_ratiosim=mean(m$ratiodiff),
        sd_ratiosim=sd(m$ratiodiff))
    }

    totalpropsimDBunkD8DCatvsOutbywindow <- rbind(AvsC_D8DBunkDCatvsD10, AvsC_D8DBunkDCatvsW1, AvsC_D8DBunkDCatvsW6, AvsC_D8DBunkDCatvsPulicaria, AvsC_D8DBunkDCatvsObtusa)

    save(totalpropsimDBunkD8DCatvsOutbywindow, file="totalpropsimDBunkD8DCatvsOutbywindow_20200824.Rdata")


    #!/usr/bin/env Rscript

        ### libraries
          library(data.table)
          library(foreach)
          library(ggplot2)
          library(tidyverse)
          library(doMC)
          registerDoMC(20)

  load("totalpropsimACotherbywindow_20200824.Rdata")
  load("totalpropsimACACoutgroupbywindow_20200824.Rdata")
  load("totalpropsimDBunkD8DCatvsOutbywindow_20200824.Rdata")

  totalttest <- rbind(totalpropsimACotherbywindow, totalpropsimACACoutgroupbywindow, totalpropsimDBunkD8DCatvsOutbywindow)

  totalttest$comppond <- factor(totalttest$comppond, levels=c("AvsC_D8D8", "AvsC_DBunkDBunk", "AvsC_DCatDCat",
    "AvsC_D8DBunk", "AvsC_D8DCat", "AvsC_DBunkDCat", "AvsC_A_D10", "AvsC_C_D10", "AvsC_D8DBunkDCatvsD10",
    "AvsC_A_W1", "AvsC_C_W1", "AvsC_D8DBunkDCatvsW1", "AvsC_A_W6", "AvsC_C_W6", "AvsC_D8DBunkDCatvsW6",
    "AvsC_A_pulicaria", "AvsC_C_pulicaria", "AvsC_D8DBunkDCatvsPulicaria", "AvsC_A_obtusa", "AvsC_C_obtusa",
    "AvsC_D8DBunkDCatvsObtusa"))

  totalttest$outgroup <- ifelse(totalttest$comppond=="AvsC_A_obtusa" | totalttest$comppond=="AvsC_A_pulicaria" |
  totalttest$comppond=="AvsC_C_obtusa" | totalttest$comppond=="AvsC_C_pulicaria" | totalttest$comppond=="AvsC_D8DBunkDCatvsPulicaria" |
  totalttest$comppond=="AvsC_D8DBunkDCatvsObtusa", "Between_Species", "Within_Species")

  totalttest$outgroup <- factor(totalttest$outgroup, levels=c("Within_Species", "Between_Species"))

  load("windowinfo_20200824.Rdata")
  setkey(windowinfo, window)
  setkey(totalttest, window)
  mtotalttest <- merge(windowinfo, totalttest)


  ggplot(data=totalttest[outgroup=="Within_Species"], aes(x=window, y=mean_ratiosim, color=comppond)) + geom_line() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept=1,color="red") + facet_wrap(~comppond)


  ggplot(data=mtotalttest[comppond=="AvsC_D8D8" | comppond=="AvsC_DBunkDBunk" | comppond=="AvsC_DCatDCat" |
  comppond=="AvsC_D8DBunk" | comppond=="AvsC_D8DCat" | comppond=="AvsC_DBunkDCat"],
  aes(x=window, y=mean_ratiosim, color=chr)) + geom_line() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept=1,color="red") + facet_wrap(~comppond)


  ggplot(data=totalttest[comppond!="AvsC_A_obtusa" & comppond!="AvsC_A_pulicaria" &
  comppond!="AvsC_C_obtusa" & comppond!="AvsC_C_pulicaria" & comppond!="AvsC_D8DBunkDCatvsPulicaria" &
  comppond!="AvsC_D8DBunkDCatvsObtusa"], aes(x=comppond, y=mean_ratiosim, group=comppond)) + geom_boxplot() +
  geom_hline(yintercept=1,color="red") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
