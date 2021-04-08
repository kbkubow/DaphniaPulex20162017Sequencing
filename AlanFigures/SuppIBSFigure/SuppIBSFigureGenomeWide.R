module load  gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3
R

#!/usr/bin/env Rscript

    ### libraries
      library(data.table)
      library(foreach)
      library(ggplot2)
      library(tidyverse)
      library(doMC)
      registerDoMC(20)

### load files

  load("totalgenomewidenoognoACsub.Rdata")

### Graph
  IBSplot <- ggplot(data=totalgenomewidenoognoACsub, aes(x=type, y=mean_ratiosim, group=type)) + geom_boxplot() +
    geom_hline(yintercept=1, color="red") + theme_bw() +
    ylab("AvsC IBS/Comparison IBS") + xlab("Comparison")

  IBSplotB <- IBSplot + scale_x_discrete(labels=c("Within D8", "Within DCat/DBunk", "Between DCat/D8/DBunk", "Kilwood vs D10", "Kilwood vs W1", "Kilwood vs W6"))

  IBSplotB + coord_flip()

### How many points are above 1
  totalgenomewidenoognoACsub[type=="AvsC_D8D8" & mean_ratiosim>=1]
  #14
  totalgenomewidenoognoACsub[type=="AvsC_D8D8" & mean_ratiosim<1]
  #1756
  14/(14+1756)
  #0.007909605

  totalgenomewidenoognoACsub[type=="AvsC_WithinPond" & mean_ratiosim>=1]
  #57
  totalgenomewidenoognoACsub[type=="AvsC_WithinPond" & mean_ratiosim<1]
  #4236
  57/(57+4236)
  #0.01327743

  totalgenomewidenoognoACsub[type=="AvsC_BetweenPond" & mean_ratiosim>=1]
  #488
  totalgenomewidenoognoACsub[type=="AvsC_BetweenPond" & mean_ratiosim<1]
  #5995
  488/(488+5995)
  #0.07527379

  (14+57+488)/(14+57+488+1756+4236+5995)
  #0.04455603

  totalgenomewidenoognoACsub[type=="AvsC_D8DBunkDCatvsD10" & mean_ratiosim>=1]
  #462
  totalgenomewidenoognoACsub[type=="AvsC_D8DBunkDCatvsD10" & mean_ratiosim<1]
  #0

  totalgenomewidenoognoACsub[type=="AvsC_D8DBunkDCatvsW1" & mean_ratiosim>=1]
  #154
  totalgenomewidenoognoACsub[type=="AvsC_D8DBunkDCatvsW1" & mean_ratiosim<1]
  #0

  totalgenomewidenoognoACsub[type=="AvsC_D8DBunkDCatvsW6" & mean_ratiosim>=1]
  #308
  totalgenomewidenoognoACsub[type=="AvsC_D8DBunkDCatvsW6" & mean_ratiosim<1]
  #0
