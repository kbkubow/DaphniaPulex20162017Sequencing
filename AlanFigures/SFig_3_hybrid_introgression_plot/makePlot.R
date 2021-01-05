
### libraries
  library(LEA)
  library(foreach)
  library(data.table)
  library(ggplot2)
  library(cowplot); theme_set(theme_cowplot())
  library(patchwork)

### load datda
  load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/SFig_3_hybrid_introgression_plot/snmf_out.Rdata")
  load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/SFig_3_hybrid_introgression_plot/fstat.Rdata")

### process LEA data
  ### source functions
    source("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/LEA/plotLeaFuncs.R")

  ### based on average
    ce.dt.ag <- ce.dt[,list(mu=median(ce), sd=sd(ce), n=length(ce)), list(k)]
    ce.dt.ag[which.min(mu)]
    ce.dt[k==8][which.min(ce)]
    k.i<-8; run.i<-30

  ### absolute lowest
    ce.dt[which.min(ce)]
    k.i<-11; run.i<-21

  ### process data
    kr <- processQ(k=k.i, run=run.i, samp=samps$clone, orderby="set2", n=3)

  ### make plot objects
    q.plot <- plotQ(kr)
    out.plot <- plotQ.sp(kr)

### Process Dsuite output
  setnames(f, "p-value", "p")
  setnames(f, "f4-ratio", "f4")

  f[data=="orig", pa:=p.adjust(p)]
  f[data=="orig"][P3=="pulicaria"][order(z)][,c('P1', 'P2', 'P3', 'Dstatistic', 'pa', 'f4', 'ABBA', 'BABA'), with=F]
  f[data=="orig"][f4>0.2]

### plot
  # panel A: SNMF individual ancestry plot (structure like plot)
    snmf.ancestry <-
      ggplot(kr, aes(x=ord, y=value, fill=as.factor(variable))) +
      geom_bar(position="stack", stat="identity", width = 1) +
      facet_grid(~set, scale="free_x", space="free") +
      geom_text(data=k.i[SC=="A"], aes(x=ord, y=1.1, label="A"), size=5) +
      geom_text(data=k.i[SC=="C"], aes(x=ord, y=1.1, label="C"), size=5) +

      theme(strip.text.x = element_text(angle=90),
            panel.spacing = unit(.2, "lines"),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_x_continuous(expand = c(0,0)) +
      xlab("") +
      guides(fill=FALSE)

  ### panels B & C: Coefficients for Dpuli & Dobt ancestry
    kr[SC=="A",col:="A"]
    kr[SC=="C",col:="C"]
    setorderv(kr, "col")

    puli.k <- k.i[Species=="pulicaria"][which.max(value)]$variable
    obt.k <- k.i[Species=="obtusa"][which.max(value)]$variable

    puli.plot <-
      ggplot(data=k.i[variable==puli.k][!is.na(set)], aes(x=set, y=log10(value), color=col)) +
      geom_jitter(width = 0.25) +
      coord_flip() +
      xlab("% D. pulicaria ancestry") +
      theme(legend.position="none")

    obt.plot <-
      ggplot(data=k.i[variable==obt.k][!is.na(set)], aes(x=set, y=log10(value), color=col)) +
      geom_jitter(width = 0.25) +
      coord_flip() +
      xlab("% D. obtusa ancestry") +
      theme(legend.position="none")

  ### panel D: optimal K


  ### panel D: Dsutie analysis

  f4 <-
  ggplot() +
  geom_point(data=f[data=="orig"][P3=="pulicaria"],
              aes(x=Dstatistic, y=log10(f4),
                 color=as.factor(pa<.05)))





    ce.plot <-   ggplot() +
                 geom_line(data=ce.dt, aes(x=k, y=ce, group=run), color="grey", alpha=.75) +
                 geom_line(data=ce.dt.ag, aes(x=k, y=mu), size=1) +
                 geom_point(data=ce.dt[k==ce.dt.ag[which.min(mu)]$k][order(ce)][1], aes(x=k, y=ce), size=3, color="black") +
                 theme(legend.position="none")

    output.plot <- q.plot / ( out.plot | ce.plot)

    ggsave(output.plot, file="~/LEA_output.new.pdf")












layout <- "
AAAAAAAA
##BBBBCC
"




  +
  plot_layout(design = layout)


plot.q1 <- plotQ(processQ(k=1, samp=samps$clone, orderby="set"))
plot.q2 <- plotQ(processQ(k=2, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q3 <- plotQ(processQ(k=3, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q4 <- plotQ(processQ(k=4, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q5 <- plotQ(processQ(k=5, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q6 <- plotQ(processQ(k=6, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q7 <- plotQ(processQ(k=7, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())

plot.q1 / plot.q2 / plot.q3 / plot.q4 / plot.q5 / plot.q6 / plot.q7



hist(k[variable.x=="V6"]$value.x)
k[variable.x=="V6"][value.x>.2]


hist(log10(k[variable.x=="V3"]$value.x), breaks=1000)
hist(log10(k[variable.x=="V6"]$value.x), breaks=1000)
