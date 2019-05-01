Looking at total ephippia production for Spring 2018 individuals established in the lab.
```
  #Load libraries
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(tidyverse)

#Read in file
        epp <- fread("EpphippiaCountsfromIndividualsKeptinLab.csv")
        colnames(epp) <- c("Pond", "Clone", "LastMinusFirstObs", "TotalEp", "EpperLmFD")
        
#Look at distribution of ephippia production
        ggplot(data=epp[LastMinusFirstObs>125], aes(x=EpperLmFD)) + geom_histogram()
        ggplot(data=epp[LastMinusFirstObs>125 & Pond=="D8"], aes(x=EpperLmFD)) + geom_histogram()
        ggplot(data=epp[LastMinusFirstObs>125 & Pond=="D8"], aes(x=TotalEp)) + geom_histogram()
        ggplot(data=epp[LastMinusFirstObs>125 & Pond!="DOily" & Pond!="Doily" & Pond!="Dmud" &
          Pond!="Doak"], aes(x=EpperLmFD, fill=Pond)) + geom_histogram() + facet_wrap(~Pond)
