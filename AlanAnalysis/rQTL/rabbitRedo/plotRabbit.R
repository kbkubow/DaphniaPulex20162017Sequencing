### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)

### load posterior probabilities
  pp <- fread("data", skip="magicReconstruct-Summary,genoprob,Conditonal genotype probability", header=T, fill=T)

  setnames(pp, names(pp), as.character(pp[1,]))
  pp <- pp[grepl("genotype", marker)]

  ppl <- melt(pp[-(1:3)], id.vars="marker")
  ppl[,value2:=as.numeric(paste(value, "0", sep=""))]

  ppl[,genotype:=last(tstrsplit(marker, "_"))]
  ppl[,clone:=gsub("_genotype[0-9]{1,}", "", marker)]

### plot
  ggplot(data=ppl, aes(x=variable, y=genotype, fill=value2)) +
  geom_tile() +
  facet_wrap(~clone) +
  
