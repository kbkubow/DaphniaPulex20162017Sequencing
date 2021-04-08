module load  gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3

R

#!/usr/bin/env Rscript

### Load libraries
  library(cn.mops)
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(tidyverse)

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
  #scaffolds <- c("Scaffold_1863_HRSCAF_2081", "Scaffold_1931_HRSCAF_2197",
    "Scaffold_2158_HRSCAF_2565", "Scaffold_2158_HRSCAF_2565",
    "Scaffold_2373_HRSCAF_2879", "Scaffold_6786_HRSCAF_7541",
    "Scaffold_7757_HRSCAF_8726", "Scaffold_9197_HRSCAF_10753",
    "Scaffold_9198_HRSCAF_10754", "Scaffold_9199_HRSCAF_10755",
    "Scaffold_9200_HRSCAF_10757", "Scaffold_9201_HRSCAF_10758")

  bamDataRangesF <- getReadCountsFromBAM(BAMFiles, refSeqName=varA, WL=750)

  resF <- cn.mops(bamDataRangesF)

  save(bamDataRangesF, file=paste(varA, "_bamDataRangesF_20200406.Rdata", sep=""))

  save(resF, file=paste(varA, "altthresh_resF_20200406.Rdata", sep=""))

### Load in data files

  cnmopsfiles <- list.files(pattern="resF_20200406.Rdata$")

  cnmopstotalregions <- foreach(i=1:length(cnmopsfiles), .combine="rbind")%do%{
    f <- cnmopsfiles[[i]]
    load(f)
    resFdt <- as.data.table(cnvr(resF))
    resFdt

  }

### Make long format

  clones <- colnames(cnmopstotalregions[,6:114])

  cnmopsregionslong <- melt(cnmopstotalregions, measure.vars=clones, variable.name="clone", value.name="ploidy")

  cnmopsregionslong$clone <- str_replace(cnmopsregionslong$clone, "_finalmap_mdup.bam", "")

  save(cnmopsregionslong, file="cnmopsregionslong.Rdata")

  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  sc$population <- str_replace(sc$population, "Dcat", "DCat")
  sc$clone <- str_replace(sc$clone, "Dcat", "DCat")

  setkey(cnmopsregionslong, clone)
  setkey(sc, clone)
  mcnmopsregionslong <- merge(cnmopsregionslong, sc)

  mcnmopsregionslong$type <- ifelse(mcnmopsregionslong$AxCF1Hybrid==1, "AxCF1", ifelse(
    mcnmopsregionslong$SC=="selfedC", "selfedC", ifelse(mcnmopsregionslong$SC=="A", "A", ifelse(
    mcnmopsregionslong$SC=="C", "C", "selfedA"
    ))
  ))

  ggplot(data=mcnmopsregionslong[seqnames=="Scaffold_2158_HRSCAF_2565" &
    start==627751 & end==633000], aes(x=type, y=ploidy)) + geom_count()
