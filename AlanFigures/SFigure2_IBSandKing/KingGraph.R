#!/usr/bin/env Rscript

### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(tidyverse)
  library(patchwork)

### Load file
  load("kinshipscm2.Rdata")
  load("ibs.longunique.Rdata")

### Graph
King <- ggplot(data=kinshipscm2[type=="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship)) + geom_point() +
  geom_point(data=kinshipscm2[type!="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship,
  color=factor(type, labels=c("Clonal", "Selfed parent-offspring",
  "Outcrossed parent-offspring", "Full-siblings", "A vs C")))) + labs(color="Type") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

IBS <- ggplot(data=ibs.longunique, aes(x=IBS)) + geom_histogram(binwidth=0.001) +
  geom_vline(xintercept = 0.965, color="red") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

IBSplusKing <- IBS + King
