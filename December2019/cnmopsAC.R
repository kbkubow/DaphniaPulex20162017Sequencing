module load  gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3

R

#!/usr/bin/env Rscript

### Load libraries
  library(cn.mops)
  library(data.table)

  ############
  ### args ###
  ############

  args=commandArgs(trailingOnly=TRUE)
  varA=args[1]

   varA

### Load SC file
  #sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  #scsub <- sc[LabGenerated==0 & Nonindependent==0 & Species=="pulex"]
  #scsubb <- scsub[SC=="A" | SC=="C" | AxCF1Hybrid==1]

### Load bam files
  BAMFiles <- list.files(pattern=".bam$")
  #BAMFiles <- list.files(path="/project/berglandlab/Karen/MappingDec2019/bams/PulexBams",
    pattern="bam$", full.names=TRUE)
  #BAMfilesB <- paste(
    scsubb$clone, "_finalmap_mdup.bam", sep="")
  #BAMfilesC <- BAMfilesB[1:5]
  #BAMfilesB <- str_replace(BAMfilesB, "May_2017_DBunk_549B", "May_2017_DBunk_549")
  bamDataRangesA <- getReadCountsFromBAM(BAMFiles, refSeqName=varA)
  bamDataRangesB <- getReadCountsFromBAM(BAMFiles, refSeqName=varA, WL=5000)
  bamDataRangesC <- getReadCountsFromBAM(BAMFiles, refSeqName=varA, WL=2500)
  bamDataRangesD <- getReadCountsFromBAM(BAMFiles, refSeqName=varA, WL=10000)

  resA <- cn.mops(bamDataRangesA)
  resB <- cn.mops(bamDataRangesB)
  resC <- cn.mops(bamDataRangesC)
  resD <- cn.mops(bamDataRangesD)

  resA <- calcIntegerCopyNumbers(resA)
  resB <- calcIntegerCopyNumbers(resB)
  resC <- calcIntegerCopyNumbers(resC)
  resD <- calcIntegerCopyNumbers(resD)

  #CNVSa <- as.data.table(cnvs(resA))
  #CNVSb <- as.data.table(cnvs(resB))
  #CNVSc <- as.data.table(cnvs(resC))
  #CNVSd <- as.data.table(cnvs(resD))

  save(bamDataRangesA, file=paste(varA, "_bamDataRangesA.Rdata")
  save(bamDataRangesB, file=paste(varA, "_bamDataRangesB.Rdata")
  save(bamDataRangesC, file=paste(varA, "_bamDataRangesC.Rdata")
  save(bamDataRangesD, file=paste(varA, "_bamDataRangesD.Rdata")

  save(resA, file=paste(varA, "_resA.Rdata"))
  save(resB, file=paste(varA, "_resB.Rdata"))
  save(resC, file=paste(varA, "_resC.Rdata"))
  save(resD, file=paste(varA, "_resD.Rdata"))
