#!/usr/bin/env Rscript

### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(tidyverse)

  ############
  ### args ###
  ############

  args=commandArgs(trailingOnly=TRUE)
  varA=args[1]

### Let's work through one input file first
  load(paste("m_IBSbyslidingwindow_250000_10000_withpulicariaandobtusa_20200629_", varA, ".Rdata", sep=""))

# Load superclone file
  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  sc$population <- str_replace(sc$population, "Dcat", "DCat")

# Add in superclone information
  scsubA <- data.table(cloneA=sc$clone, SCA=sc$SC, popA=sc$population, yearA=sc$year, medrdA=sc$medrd, speciesA=sc$Species, sexA=sc$Sex)
  scsubB <- data.table(cloneB=sc$clone, SCB=sc$SC, popB=sc$population, yearB=sc$year, medrdB=sc$medrd, speciesB=sc$Species, sexB=sc$Sex)
  setkey(scsubA, cloneA)
  setkey(scsubB, cloneB)
  setkey(m, cloneB)
  mtmp <- merge(m, scsubB)
  setkey(mtmp, cloneA)
  msc <- merge(mtmp, scsubA)

# Now add in new composite superclone info columns, start with making new OO SCs
  msc$SCA_B <- ifelse(msc$SCA=="OO", paste(msc$SCA, msc$cloneA, sep="_"), msc$SCA)
  msc$SCB_B <- ifelse(msc$SCB=="OO", paste(msc$SCB, msc$cloneB, sep="_"), msc$SCB)
  msc$compinfo <- paste(msc$SCA_B, msc$popA, msc$SCB_B, msc$popB, sep="_")
  msc$SCcomp <- paste(msc$SCA, msc$SCB, sep="_")

  # And I will make a file to compare all non A and C lineages from the three focal ponds to the outgroup species

      mscAC <- msc[SCA=="A" | SCA=="C" | SCB=="A" | SCB=="C"]
      mscACnoObtusaPulicaria <- mscAC[SCA!="G" & SCA!="P" & SCB!="G" & SCB!="P"]
      mscACnoObtusaPulicariaotheronly <- mscACnoObtusaPulicaria[popA!="D10" & popB!="D10" &
        popA!="DLily" & popB!="DLily" & popA!="DOil" & popB!="DOil" & popA!="Dramp" &
        popB!="Dramp" & popA!="W1" & popB!="W1" & popA!="W6" & popB!="W6" & popA!="DMud" & popB!="DMud"]
      mscACnoObtusaPulicariaotheronly <- mscACnoObtusaPulicariaotheronly[SCcomp!="A_C" & SCcomp!="C_A"]
      mscACnoObtusaPulicariaotheronly$ACvspond <- ifelse(mscACnoObtusaPulicariaotheronly$SCA=="A" |
        mscACnoObtusaPulicariaotheronly$SCB=="A", "D8DBunkDCatvsA", "D8DBunkDCatvsC")
      mscACnoObtusaPulicariaotheronly$SCA_C <- ifelse(
        mscACnoObtusaPulicariaotheronly$ACvspond=="D8DBunkDCatvsA", "A", "C")
      mscACnoObtusaPulicariaotheronly$SCB_C <- ifelse(
        mscACnoObtusaPulicariaotheronly$SCB_B=="A" | mscACnoObtusaPulicariaotheronly$SCB_B=="C",
        mscACnoObtusaPulicariaotheronly$SCA_B, mscACnoObtusaPulicariaotheronly$SCB_B)

      D8DBunkDCatvsAC.ag <- mscACnoObtusaPulicariaotheronly[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
        list(window, chr, start, stop, speciesA, speciesB, SCA_C, SCB_C, ACvspond)]

      save(D8DBunkDCatvsAC.ag, file=paste("D8DBunkDCatvsAC.ag_", varA, ".Rdata", sep=""))

      inputfilesOtherAC <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="D8DBunkDCatvsAC.ag_")

      totalibsOtherAC <- foreach(i=1:length(inputfilesOtherAC), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesOtherAC[i], sep="")
        load(f)
        D8DBunkDCatvsAC.ag
      }

      totalibsOtherAC_A <- totalibsOtherAC[ACvspond=="D8DBunkDCatvsA"]
      totalibsOtherAC_C <- totalibsOtherAC[ACvspond=="D8DBunkDCatvsC"]
      totalibsOtherAC_Asub <- totalibsOtherAC_A[, c("window", "chr", "start", "stop", "SCB_C", "meanIBS")]
      colnames(totalibsOtherAC_Asub) <- c("window", "chr", "start", "stop", "SCB_C", "meanAIBS")
      totalibsOtherAC_Csub <- totalibsOtherAC_C[, c("window", "chr", "start", "stop", "SCB_C", "meanIBS")]
      colnames(totalibsOtherAC_Csub) <- c("window", "chr", "start", "stop", "SCB_C", "meanCIBS")
      setkey(totalibsOtherAC_Asub, window, chr, start, stop, SCB_C)
      setkey(totalibsOtherAC_Csub, window, chr, start, stop, SCB_C)
      m <- merge(totalibsOtherAC_Asub, totalibsOtherAC_Csub)

      totest <- unique(m$SCB_C)

      AvsC_D8D8 <- foreach(i=1:length(totest), .combine="rbind")%do%{
        print(paste("Getting data: ", i, sep=""))
        SCB <- totest[i]
        tmp <- m[SCB_C==SCB]
        ttest <- t.test(tmp$meanAIBS, tmp$meanCIBS, paired=TRUE)
        stats <- data.table(compsc=SCB, t=ttest$statistic, df=ttest$parameter,
        p=ttest$p.value, CIlow=ttest$conf.int[1], CIhigh=ttest$conf.int[2],
        meandiff=ttest$estimate[1])
      }

      AvsC_D8D8$sig0.5 <- ifelse(AvsC_D8D8$p < 0.05, 1, 0)

      ggplot(data=AvsC_D8D8, aes(x=t, fill=as.factor(sig0.5))) + geom_histogram(binwidth=5)
      ggplot(data=AvsC_D8D8, aes(x=meandiff, fill=as.factor(sig0.5))) + geom_histogram(binwidth=5)
