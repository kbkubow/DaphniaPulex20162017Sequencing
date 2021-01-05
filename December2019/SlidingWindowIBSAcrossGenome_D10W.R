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

      inputfilesOtherAC <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="D8DBunkDCatvsAC.ag_")

      totalibsOtherAC <- foreach(i=1:length(inputfilesOtherAC), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesOtherAC[i], sep="")
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
      vsD10W <- total[type=="vsD10" | type=="vsW1" | type=="vsW6"]
      window <- unique(vsD10W$window)
      ACvspondu <- unique(vsD10W$ACvspondB)

      focalbywindowbrokendownD10W1 <- foreach(n=1:length(ACvspondu), .combine="rbind")%do%{
          AC <- ACvspondu[n]
          subsetfocal <- vsD10W[ACvspondB==AC]
          focalbywindow <- foreach(i=1:length(window), .combine="rbind")%do%{
            w <- window[i]
            tmp <- subsetfocal[window==w]
            q <- as.data.table(quantile(tmp$meanIBS,probs=c(0.025,0.975)))
            tmp2 <- data.table(window=w, ACvspond=unique(tmp$ACvspondB), low95=q[1], high95=q[2], AvsCIBS=AvsC$meanIBS[window==w])
            tmp2
          }
          focalbywindow
      }

      save(focalbywindowbrokendownD10W1, file="focalbywindowbrokendownD10W1_20201015.Rdata")




  ### libraries
        library(data.table)
        library(foreach)
        library(ggplot2)
        library(tidyverse)

  ### Load IBS output file, window file, and merge

      load("focalbywindow_20201008.Rdata")
      load("windowinfo_20200824.Rdata")
      PA42 <- fread("5kbchrassignHiCnew.csv")

      setkey(focalbywindow, window)
      setkey(windowinfo, window)
      mfocalbywindow <- merge(focalbywindow, windowinfo)

      PA42main <- PA42[MainChromosome=="Y"]
      colnames(PA42main) <- c("chr", "Length", "PA42chr", "MainChr")
      subPA42main <- PA42main[, c("chr", "PA42chr"), with=FALSE]

      setkey(subPA42main, chr)
      setkey(mfocalbywindow, chr)
      m2focalbywindow <- merge(mfocalbywindow, subPA42main)

      m2focalbywindow$PA42chr <- factor(m2focalbywindow$PA42chr, levels=c("1", "2", "3", "4",
        "5", "6", "7", "8", "9", "10", "11", "12"))

      m2focalbywindow$window <- as.factor(m2focalbywindow$window)

      m2focalbywindow$tograph <- c("all")

  ### Graph IBS
      ggplot() + geom_line(data=m2focalbywindow, aes(x=reorder(window, PA42chr), y=AvsCIBS, color=as.factor(PA42chr), group=tograph)) +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      ggplot() + geom_line(data=m2focalbywindow, aes(x=reorder(window, PA42chr), y=AvsCIBS, color=as.factor(PA42chr), group=tograph), size=1) +
      geom_line(data=m2focalbywindow, aes(x=reorder(window, PA42chr), y=low95.V1, group=tograph), color="grey") +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      IBSwithconfidence <- ggplot() + geom_line(data=m2focalbywindow, aes(x=reorder(window, PA42chr), y=AvsCIBS, color=as.factor(PA42chr), group=tograph), size=1) +
      geom_line(data=m2focalbywindow, aes(x=reorder(window, PA42chr), y=low95.V1, group=tograph, alpha=0.75), size=0.35) +
      geom_line(data=m2focalbywindow, aes(x=reorder(window, PA42chr), y=high95.V1, group=tograph, alpha=0.75), size=0.35) +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      ggsave(IBSwithconfidence, file="IBSwithconfidence_AvsCvsfocal_20201009.pdf")
