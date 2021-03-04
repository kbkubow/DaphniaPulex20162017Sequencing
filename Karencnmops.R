#!/usr/bin/env Rscript

### Load libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(cnmops)

### Load files
  load("/Users/kbkubow/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/STAR/DESeq_star.output")
  gff <- fread("/Users/kbkubow/Box Sync/UVA/Daphnia/ReferenceGenome/Dovetail/newHiC/Daphnia.aed.0.6.gff")
  colnames(gff) <- c("chr", "maker", "type", "start", "stop", "V6", "V7", "V8", "Notes")
  tmp <- sapply(strsplit(gff$Notes, ";"), getElement, 1)
  tmp2 <- sapply(strsplit(tmp, "="), getElement, 2)
  tmp3 <- sapply(strsplit(tmp2, ":"), getElement, 1)
  gff$GeneId <- tmp3

  resEfiles <- list.files(pattern="resE.Rdata")

  totalresE <- foreach(i=1:length(resEfiles), .combine="rbind")%do%{
    f=paste("/Users/kbkubow/Box Sync/UVA/Daphnia/InitialManuscript/cnmops/", resEfiles[i], sep="")
    load(f)
    cnvr <- as.data.table(cnvr(resE))
    cnvr
  }

  resFfiles <- list.files(pattern="resF.Rdata")

  totalresF <- foreach(i=1:length(resFfiles), .combine="rbind")%do%{
    f=paste("/Users/kbkubow/Box Sync/UVA/Daphnia/InitialManuscript/cnmops/", resFfiles[i], sep="")
    load(f)
    cnvr <- as.data.table(cnvr(resF))
    cnvr
  }

  resFaltfiles <- list.files(path="./altthres/", pattern="resF.Rdata")

  totalaltresF <- foreach(i=1:length(resFaltfiles), .combine="rbind")%do%{
    f=paste("/Users/kbkubow/Box Sync/UVA/Daphnia/InitialManuscript/cnmops/altthres/", resFaltfiles[i], sep="")
    load(f)
    cnvr <- as.data.table(cnvr(resF))
    cnvr
  }

  resGfiles <- list.files(pattern="resG.Rdata")

  totalresG <- foreach(i=1:length(resFfiles), .combine="rbind")%do%{
    f=paste("/Users/kbkubow/Box Sync/UVA/Daphnia/InitialManuscript/cnmops/", resGfiles[i], sep="")
    load(f)
    cnvr <- as.data.table(cnvr(resG))
    cnvr
  }


### Pull out highly differentially expressed genes
  # Negative log2foldchange indicate higher expression in A
  res.dthigh <- res.dt[abs(log2FoldChange) > 3]

  setkey(res.dthigh, GeneId)
  setkey(gff, GeneId)
  m <- merge(res.dthigh, gff)

  setkey(res.dt, GeneId)
  setkey(gff, GeneId)
  m2 <- merge(res.dt, gff)


  highlyexpinCNV <- foreach(i=1:dim(totalresF)[1], .combine="rbind")%do%{
    scaff <- totalresF$seqnames[[i]]
    st <- totalresF$start[[i]] - 1000
    ed <- totalresF$end[[i]] + 1000
    D8179p <- totalresF$April_2017_D8_179_finalmap_mdup.bam[[i]]
    D8222p <- totalresF$April_2017_D8_222_finalmap_mdup.bam[[i]]
    D8349p <- totalresF$April_2017_D8_349_finalmap_mdup.bam[[i]]
    D8515p <- totalresF$May_2017_D8_515_finalmap_mdup.bam[[i]]
    tmp <- m[chr==scaff & start > st & start < ed | chr==scaff & stop > st & stop < ed]
    #tmp <- m[chr==scaff & start > st & stop < ed]
    tmp$D8179 <- D8179p
    tmp$D8222 <- D8222p
    tmp$D8349 <- D8349p
    tmp$D8515 <- D8515p
    tmp
  }


  expinCNV <- foreach(i=1:dim(m2)[1], .combine="rbind")%do%{
    scaff <- m2$chr[[i]]
    st <- m2$start[[i]]
    ed <- m2$stop[[i]]
    tmpres <- totalresF[seqnames==scaff & st > start-1000 & st < end+1000 | seqnames==scaff & ed > start-1000 & ed < end+1000]
    D8179p <- mean(as.numeric(str_replace(tmpres$April_2017_D8_179_finalmap_mdup.bam, "CN", "")))
    D8222p <- mean(as.numeric(str_replace(tmpres$April_2017_D8_222_finalmap_mdup.bam, "CN", "")))
    D8349p <- mean(as.numeric(str_replace(tmpres$April_2017_D8_349_finalmap_mdup.bam, "CN", "")))
    D8515p <- mean(as.numeric(str_replace(tmpres$May_2017_D8_515_finalmap_mdup.bam, "CN", "")))
    tmp <- m2[i]
    tmp$D8179 <- D8179p
    tmp$D8222 <- D8222p
    tmp$D8349 <- D8349p
    tmp$D8515 <- D8515p
    tmp
  }

  expinCNV[is.na(D8179),D8179:=2]
  expinCNV[is.na(D8222),D8222:=2]
  expinCNV[is.na(D8349),D8349:=2]
  expinCNV[is.na(D8515),D8515:=2]

  expinCNV$meanCminusA <- ((expinCNV$D8515+expinCNV$D8222)/2) - ((expinCNV$D8179+expinCNV$D8349)/2)

  ggplot(data=expinCNV, aes(x=meanCminusA, y=log2FoldChange)) + geom_point()




  genesCNVs <- foreach(i=1:dim(totalresF)[1], .combine="rbind")%do%{
    scaff <- totalresF$seqnames[[i]]
    st <- totalresF$start[[i]] - 1000
    ed <- totalresF$end[[i]] + 1000
    D8179p <- totalresF$April_2017_D8_179_finalmap_mdup.bam[[i]]
    D8222p <- totalresF$April_2017_D8_222_finalmap_mdup.bam[[i]]
    D8349p <- totalresF$April_2017_D8_349_finalmap_mdup.bam[[i]]
    D8515p <- totalresF$May_2017_D8_515_finalmap_mdup.bam[[i]]
    tmp <- m2[chr==scaff & start > st & start < ed | chr==scaff & stop > st & stop < ed]
    #tmp <- m[chr==scaff & start > st & stop < ed]
    tmp$D8179 <- D8179p
    tmp$D8222 <- D8222p
    tmp$D8349 <- D8349p
    tmp$D8515 <- D8515p
    tmp
  }


  expinCNV$D8179 <- str_replace(expinCNV$D8179, "CN0", "0")
  expinCNV$D8179 <- str_replace(expinCNV$D8179, "CN1", "1")
  expinCNV$D8179 <- str_replace(expinCNV$D8179, "CN2", "2")
  expinCNV$D8179 <- str_replace(expinCNV$D8179, "CN3", "3")
  expinCNV$D8179 <- str_replace(expinCNV$D8179, "CN4", "4")
  expinCNV$D8222 <- str_replace(expinCNV$D8222, "CN0", "0")
  expinCNV$D8222 <- str_replace(expinCNV$D8222, "CN1", "1")
  expinCNV$D8222 <- str_replace(expinCNV$D8222, "CN2", "2")
  expinCNV$D8222 <- str_replace(expinCNV$D8222, "CN3", "3")
  expinCNV$D8349 <- str_replace(expinCNV$D8349, "CN0", "0")
  expinCNV$D8349 <- str_replace(expinCNV$D8349, "CN1", "1")
  expinCNV$D8349 <- str_replace(expinCNV$D8349, "CN2", "2")
  expinCNV$D8349 <- str_replace(expinCNV$D8349, "CN3", "3")
  expinCNV$D8349 <- str_replace(expinCNV$D8349, "CN4", "4")
  expinCNV$D8349 <- str_replace(expinCNV$D8349, "CN5", "5")
  expinCNV$D8349 <- str_replace(expinCNV$D8349, "CN6", "6")
  expinCNV$D8349 <- str_replace(expinCNV$D8349, "CN7", "7")
  expinCNV$D8349 <- str_replace(expinCNV$D8349, "CN8", "8")
  expinCNV$D8515 <- str_replace(expinCNV$D8515, "CN0", "0")
  expinCNV$D8515 <- str_replace(expinCNV$D8515, "CN1", "1")
  expinCNV$D8515 <- str_replace(expinCNV$D8515, "CN2", "2")
  expinCNV$D8515 <- str_replace(expinCNV$D8515, "CN3", "3")
  expinCNV$D8179 <- as.numeric(expinCNV$D8179)
  expinCNV$D8222 <- as.numeric(expinCNV$D8222)
  expinCNV$D8349 <- as.numeric(expinCNV$D8349)
  expinCNV$D8515 <- as.numeric(expinCNV$D8515)

  expinCNV$meanCminusA <- ((expinCNV$D8515+expinCNV$D8222)/2) - ((expinCNV$D8179+expinCNV$D8349)/2)

  ggplot(data=expinCNV, aes(x=meanCminusA, y=log2FoldChange)) + geom_point()

  expinCNVgenes <- expinCNV$GeneId
  totalgenes <- m2$GeneId
  genesnotinCNV <- setdiff(totalgenes, expinCNVgenes)


  5 - 73; 6 - 55;
