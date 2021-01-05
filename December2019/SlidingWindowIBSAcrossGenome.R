#!/usr/bin/env Rscript

    ### libraries
      library(data.table)
      library(foreach)
      library(ggplot2)
      library(tidyverse)

  ### Load processed IBS by window files
      inputfiles <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", pattern="mscACvsrest.ag")

      totalibsag <- foreach(i=1:length(inputfiles), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", inputfiles[i], sep="")
        load(f)
        mscACvsrest.ag
      }

      inputfilesPO <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", pattern="mscACvsrest_POBC")

      totalibsagPO <- foreach(i=1:length(inputfilesPO), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", inputfilesPO[i], sep="")
        load(f)
        mscACvsrest_POBC
      }

      totalibsagPO$ACvspond <- ifelse(totalibsagPO$ACvspond=="A_obtusa" & totalibsagPO$SCA_C=="C", "C_obtusa", ifelse(
        totalibsagPO$ACvspond=="A_pulicaria" & totalibsagPO$SCA_C=="C", "C_pulicaria", totalibsagPO$ACvspond
      ))

      inputfilesDBunkD8 <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", pattern="D8DBunkDCat.ag")

      totalibsagD8DBunk <- foreach(i=1:length(inputfilesDBunkD8), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", inputfilesDBunkD8[i], sep="")
        load(f)
        D8DBunkDCat.ag
      }


      inputfilesObtusaPulicaria <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", pattern="D8DBunkDCatvsPO.ag")

      totalibsObPul <- foreach(i=1:length(inputfilesObtusaPulicaria), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", inputfilesObtusaPulicaria[i], sep="")
        load(f)
        D8DBunkDCatvsPO.ag
      }


      inputfilesD10W <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", pattern="D8DBunkDCatvsD10W.ag")

      totalibsD10W <- foreach(i=1:length(inputfilesD10W), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", inputfilesD10W[i], sep="")
        load(f)
        D8DBunkDCatvsD10W.ag
      }

      inputfilesOtherAC <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", pattern="D8DBunkDCatvsAC.ag_")

      totalibsOtherAC <- foreach(i=1:length(inputfilesOtherAC), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/slidingwindowIBS/", inputfilesOtherAC[i], sep="")
        load(f)
        D8DBunkDCatvsAC.ag
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

      #load("D8DBunkDCatvsAC.ag_G.Rdata")
      #totalibsOtherAC <- D8DBunkDCatvsAC.ag

      total <- rbind(totalibsag, totalibsagPO, totalibsagD8DBunk, totalibsD10W, totalibsObPul, totalibsOtherAC)

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
        total$ACvspondB=="C_D10" | total$ACvspondB=="D8DBunkDCatvsD10", "vsD10", ifelse(total$ACvspondB=="A_W1" |
        total$ACvspondB=="C_W1" | total$ACvspondB=="D8DBunkDCatvsW1", "vsW1", ifelse(total$ACvspondB=="A_W6" |
        total$ACvspondB=="C_W6" | total$ACvspondB=="D8DBunkDCatvsW6", "vsW6", ifelse(total$ACvspondB=="A_pulicaria" |
        total$ACvspondB=="C_pulicaria" | total$ACvspondB=="D8DBunkDCatvsPulicaria" , "vsPulicaria", ifelse(
        total$ACvspondB=="D8DBunkDCatvsA" | total$ACvspondB=="D8DBunkDCatvsC", "AvsCvsother", "vsObtusa"
      ))))))))

      AvsC <- total[type=="AvsC"]
      focalpond <- total[type=="betweenpond" | type=="withinpond"]
      window <- unique(focalpond$window)
      ACvspondu <- unique(focalpond$ACvspond)

      #focalbywindow <- foreach(i=1:length(window), .combine="rbind")%do%{
      #  w <- window[i]
      #  tmp <- subsetfocal[window==w]
      #  q <- as.data.table(quantile(tmp$meanIBS,probs=c(0.025,0.975)))
      #  tmp2 <- data.table(window=w, low95=q[1], high95=q[2], AvsCIBS=AvsC$meanIBS[window==w])
      #  tmp2
      #}

      save(focalbywindow, file="focalbywindow_20201008.Rdata")

      focalbywindowbrokendown <- foreach(n=1:length(ACvspondu), .combine="rbind")%do%{
          AC <- ACvspondu[n]
          subsetfocal <- focalpond[ACvspond==AC]
          focalbywindow <- foreach(i=1:length(window), .combine="rbind")%do%{
            w <- window[i]
            tmp <- subsetfocal[window==w]
            q <- as.data.table(quantile(tmp$meanIBS,probs=c(0.025,0.975)))
            tmp2 <- data.table(window=w, ACvspond=unique(tmp$ACvspond), low95=q[1], high95=q[2], AvsCIBS=AvsC$meanIBS[window==w])
            tmp2
          }
          focalbywindow
      }

      save(focalbywindowbrokendown, file="focalbywindowbrokendown_20201015.Rdata")




  ### libraries
        library(data.table)
        library(foreach)
        library(ggplot2)
        library(tidyverse)

  ### Load IBS output file, window file, and merge

      load("focalbywindow_20201008.Rdata")
      load("focalbywindowbrokendown_20201015.Rdata")
      load("focalbywindowbrokendownD10W1_20201015.Rdata")
      load("focalbywindowbrokendownObPul_20201016.Rdata")
      load("windowinfo_20200824.Rdata")
      PA42 <- fread("5kbchrassignHiCnew.csv")

      totalbywindow <- rbind(focalbywindowbrokendown, focalbywindowbrokendownD10W1)

      setkey(totalbywindow, window)
      setkey(windowinfo, window)
      mfocalbywindow <- merge(totalbywindow, windowinfo)

      setkey(focalbywindowbrokendownObPul, window)
      setkey(windowinfo, window)
      mfocalbywindowObPul <- merge(focalbywindowbrokendownObPul, windowinfo)


      PA42main <- PA42[MainChromosome=="Y"]
      colnames(PA42main) <- c("chr", "Length", "PA42chr", "MainChr")
      subPA42main <- PA42main[, c("chr", "PA42chr"), with=FALSE]

      setkey(subPA42main, chr)
      setkey(mfocalbywindow, chr)
      m2focalbywindow <- merge(mfocalbywindow, subPA42main)

      setkey(subPA42main, chr)
      setkey(mfocalbywindowObPul, chr)
      m2focalbywindowObPul <- merge(mfocalbywindowObPul, subPA42main)

      m2focalbywindow$PA42chr <- factor(m2focalbywindow$PA42chr, levels=c("1", "2", "3", "4",
        "5", "6", "7", "8", "9", "10", "11", "12"))

      m2focalbywindow$window <- as.factor(m2focalbywindow$window)

      m2focalbywindow$tograph <- c("all")

      m2focalbywindow$ACvspond <- factor(m2focalbywindow$ACvspond, levels=c("D8_D8", "DBunk_DBunk", "DCat_DCat",
       "D8_DBunk", "D8_DCat", "DBunk_DCat", "A_D10", "C_D10", "D8DBunkDCatvsD10", "A_W1", "C_W1",
       "D8DBunkDCatvsW1", "A_W6", "C_W6", "D8DBunkDCatvsW6"))

      m2focalbywindowObPul$PA42chr <- factor(m2focalbywindowObPul$PA42chr, levels=c("1", "2", "3", "4",
       "5", "6", "7", "8", "9", "10", "11", "12"))

      m2focalbywindowObPul$window <- as.factor(m2focalbywindowObPul$window)

      m2focalbywindowObPul$tograph <- c("all")

      m2focalbywindowObPul$ACvspond <- factor(m2focalbywindowObPul$ACvspond, levels=c("A_pulicaria",
       "C_pulicaria", "D8DBunkDCatvsPulicaria", "A_obtusa", "C_obtusa", "D8DBunkDCatvsObtusa"))

  ### Graph IBS

  IBSslidingwindowacrossgenomeD8D8only <- ggplot(m2focalbywindow[ACvspond=="D8_D8"], 
          aes(window)) +
          geom_ribbon(aes(ymin=low95.V1, ymax=high95.V1, group=tograph), fill="grey80") +
          geom_line(aes(y=AvsCIBS, group=tograph, color=as.factor(chr))) +
          theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          facet_grid(ACvspond~chr, scales="free_x")

  ggsave(IBSslidingwindowacrossgenomeD8D8only, file="IBSslidingwindowacrossgenomeD8D8only.pdf")

    IBSslidingwindowacrossgenome <- ggplot(m2focalbywindow[ACvspond=="D8_D8" | ACvspond=="DBunk_DBunk" |
              ACvspond=="DCat_DCat" | ACvspond=="D8_DBunk" | ACvspond=="D8_DCat" | ACvspond=="DBunk_DCat"],
              aes(window)) +
              geom_ribbon(aes(ymin=low95.V1, ymax=high95.V1, group=tograph), fill="grey80") +
              geom_line(aes(y=AvsCIBS, group=tograph, color=as.factor(chr))) +
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
              facet_grid(ACvspond~chr, scales="free_x")

      ggsave(IBSslidingwindowacrossgenome, file="IBSslidingwindowacrossgenome.pdf")

      IBSslidingwindowacrossgenomeD10W <- ggplot(m2focalbywindow[ACvspond=="A_D10" | ACvspond=="C_D10" |
              ACvspond=="D8DBunkDCatvsD10" | ACvspond=="A_W1" | ACvspond=="C_W1" | ACvspond=="D8DBunkDCatvsW1" |
              ACvspond=="A_W6" | ACvspond=="C_W6" | ACvspond=="D8DBunkDCatvsW6"], aes(window)) +
              geom_ribbon(aes(ymin=low95.V1, ymax=high95.V1, group=tograph), fill="grey80") +
              geom_line(aes(y=AvsCIBS, group=tograph, color=as.factor(chr))) +
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
              facet_grid(ACvspond~chr, scales="free_x")

      ggsave(IBSslidingwindowacrossgenomeD10W, file="IBSslidingwindowacrossgenomeD10W.pdf")

      IBSslidingwindowacrossgenomeD10Wwlines <- ggplot(m2focalbywindow[ACvspond=="A_D10" | ACvspond=="C_D10" |
              ACvspond=="D8DBunkDCatvsD10" | ACvspond=="A_W1" | ACvspond=="C_W1" | ACvspond=="D8DBunkDCatvsW1" |
              ACvspond=="A_W6" | ACvspond=="C_W6" | ACvspond=="D8DBunkDCatvsW6"], aes(window)) +
              geom_ribbon(aes(ymin=low95.V1, ymax=high95.V1, group=tograph), fill="grey80") +
              geom_line(aes(y=low95.V1, group=tograph), color="grey80") +
              geom_line(aes(y=high95.V1, group=tograph), color="grey80") +
              geom_line(aes(y=AvsCIBS, group=tograph, color=as.factor(chr))) +
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
              facet_grid(ACvspond~chr, scales="free_x")

      ggsave(IBSslidingwindowacrossgenomeD10Wwlines, file="IBSslidingwindowacrossgenomeD10Wwlines.pdf")

      IBSslidingwindowacrossgenomeObPul <- ggplot(m2focalbywindowObPul, aes(window)) +
              geom_ribbon(aes(ymin=low95.V1, ymax=high95.V1, group=tograph), fill="grey80") +
              geom_line(aes(y=AvsCIBS, group=tograph, color=as.factor(chr))) +
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
              facet_grid(ACvspond~chr, scales="free_x")

      ggsave(IBSslidingwindowacrossgenomeObPul, file="IBSslidingwindowacrossgenomeObPul.pdf")

      IBSslidingwindowacrossgenomeObPulwlines <- ggplot(m2focalbywindowObPul, aes(window)) +
              geom_ribbon(aes(ymin=low95.V1, ymax=high95.V1, group=tograph), fill="grey80") +
              geom_line(aes(y=low95.V1, group=tograph), color="grey80") +
              geom_line(aes(y=high95.V1, group=tograph), color="grey80") +
              geom_line(aes(y=AvsCIBS, group=tograph, color=as.factor(chr))) +
              theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
              facet_grid(ACvspond~chr, scales="free_x")

      ggsave(IBSslidingwindowacrossgenomeObPulwlines, file="IBSslidingwindowacrossgenomeObPulwlines.pdf")
