#!/usr/bin/env Rscript

### libraries
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(tidyverse)
        library(lme4)


### Load data file
  neo <- fread("FinalNeoNotAdj.csv")

  neopretrmt <- neo[, c("Clone", "MomJar", "Neonate", "Neo2ndGen", "Isolated", "Clutch1", "Male1", "Female1", "Trmt"), with=FALSE]
  neopretrmtMC <- neopretrmt[Trmt=="MF" | Trmt=="C"][Clone!="SW4" & Clone!="AxoMorse"]
  neopretrmtMC$propmale <- neopretrmtMC$Male1/(neopretrmtMC$Male1+neopretrmtMC$Female1)
  neopretrmtMC$N <- (neopretrmtMC$Male1+neopretrmtMC$Female1)
  neopretrmtMC[,PrePost:=("Pre")]
  neoposttrmt <- neo[, c("Clone", "MomJar", "Neonate", "Neo2ndGen", "Isolated", "Clutch1", "Week3PlusMales", "Week3PlusFemales", "Trmt"), with=FALSE]
  neoposttrmtMC <- neoposttrmt[Trmt=="MF" | Trmt=="C"][Clone!="SW4" & Clone!="AxoMorse"]
  neoposttrmtMC$propmale <- neoposttrmtMC$Week3PlusMales/(neoposttrmtMC$Week3PlusMales+neoposttrmtMC$Week3PlusFemales)
  neoposttrmtMC$N <- (neoposttrmtMC$Week3PlusMales+neoposttrmtMC$Week3PlusFemales)
  neoposttrmtMC[,PrePost:=("Post")]
  dat <- rbind(neopretrmtMC[,c("Clone", "MomJar", "Trmt", "PrePost", "propmale", "N"), with=F],
               neoposttrmtMC[,c("Clone", "MomJar", "Trmt", "PrePost", "propmale", "N"), with=F])
  dat[grepl("D8515|D8222", Clone), SC:="C"]
  dat[!grepl("D8515|D8222", Clone), SC:="A"]

  dat.ag <- dat[,list(propmale=sum(propmale*N, na.rm=T)/sum(N, na.rm=T), N=sum(N, na.rm=T)), list(Clone, Trmt, PrePost, SC)]
  dat.ag[,se:=sqrt((propmale*(1-propmale))/N)]
  dat.ag[,lci:=propmale-1.96*se]
  dat.ag[,uci:=propmale+1.96*se]
  dat.ag[,gr:=paste(Clone, Trmt, sep=":")]
  dat.ag$PrePost <- factor(dat.ag$PrePost, levels=c("Pre", "Post"))

  Methyl <- ggplot(data=dat.ag, aes(x=PrePost, group=gr, y=propmale, color=SC, linetype=Clone)) +
    geom_line(position=position_dodge(0.3)) +
    scale_size_manual(values = c(control=0.5,MF=1)) +
    geom_errorbar(aes(ymin=lci, ymax=uci), , linetype=1, width=.3, position=position_dodge(0.3)) +
    geom_point(data=dat.ag, aes(x=PrePost, group=gr, y=propmale, color=SC, shape=Trmt),
    position=position_dodge(0.3), size=2) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    xlab("Pre/Post Treatment") + ylab("Prop Male") +
    guides(color=FALSE) + theme(legend.position = c(0.21, 0.8)) + labs(shape="Treatment") +
    guides(linetype = FALSE)

  t1 <- glmer(propmale~1+(1|Clone)+PrePost*Trmt, data=dat, weights=N, family=binomial)
  t1 <- glmer(propmale~1+(1|Clone)+SC, data=dat[Trmt=="MF" & PrePost=="Post"], weights=N, family=binomial)
