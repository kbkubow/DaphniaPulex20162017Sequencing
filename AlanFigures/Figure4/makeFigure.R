### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(patchwork)

### setwd
  setwd("/Users/kbkubow/Documents/GitHub")

### load F1 cross data (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")

### load in PoolSeq data (made by `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLSeqR.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")

### load in F1 mapping data
  #load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")
  #load("~/lme4qtl_output.AxC.long.Rdata")
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.SUMMARIZED.Rdata")

### load in overlap tests
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/overlap_perm.Rdata")

### Load in PA42 chr key
  PA42 <- fread("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/5kbchrassignHiCnew.csv")

### male proportion plot
  male$gr <- ifelse(male$SC=="selfedA", "A", ifelse(male$SC=="B", "B", male$gr))
  male.ag <- male[!is.na(gr),list(propmale=sum(Males)/sum(NewTotal),
                                  N=sum(NewTotal)),
                    list(clone, gr)]
  male.ag[,se:=sqrt((propmale*(1-propmale))/N)]

  male.ag[,lci:=propmale-1.96*se]
  male.ag[,uci:=propmale+1.96*se]

  male.ag$gr <- factor(male.ag$gr, levels=c("A", "AxC", "C", "CxC", "B"))

  male.plot <- ggplot(data=male.ag, aes(x=gr, y=propmale, color=gr)) +
  geom_linerange(aes(ymin=lci, ymax=uci), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Male proportion") +theme_bw()

### Fill Rate
  epp$gr <- ifelse(epp$SC=="selfedA", "A", ifelse(epp$SC=="B", "B", epp$gr))
  epp$counts <- c(1)

  epp.ag <- epp[!is.na(gr),list(fillrate=sum(fill*TotalEppB, na.rm=T)/sum(TotalEppB),
                                N=sum(TotalEppB), meanEMB=mean(TotalEmb),
                                sdEMB=sd(TotalEmb), sampsize=sum(counts),
                                meanEpp=mean(TotalEpp), sdEpp=sd(TotalEpp)),
                  list(clone, gr)]

  epp.ag[,se:=sqrt((fillrate*(1-fillrate))/N)]
  epp.ag[,lci:=fillrate-1.96*se]
  epp.ag[,uci:=fillrate+1.96*se]

  epp.ag[,seemb:=sdEMB/(sqrt(sampsize))]
  epp.ag[,lciemb:=meanEMB-1.96*seemb]
  epp.ag[,uciemb:=meanEMB+1.96*seemb]

  epp.ag[,seepp:=sdEpp/(sqrt(sampsize))]
  epp.ag[,lciepp:=meanEpp-1.96*seepp]
  epp.ag[,uciepp:=meanEpp+1.96*seepp]


  epp.ag$gr <- factor(epp.ag$gr, levels=c("A", "AxC", "C", "CxC", "B"))

  epp.plot <- ggplot(data=epp.ag, aes(x=gr, y=fillrate, color=gr)) +
  geom_linerange(aes(ymin=lci, ymax=uci), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Ephippial fill rate") +
  theme_bw()

  emb.plot <- ggplot(data=epp.ag, aes(x=gr, y=meanEMB, color=gr)) +
  geom_linerange(aes(ymin=lciemb, ymax=uciemb), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Embryo production") +
  theme_bw()

  eppnum.plot <- ggplot(data=epp.ag, aes(x=gr, y=meanEpp, color=gr)) +
  geom_linerange(aes(ymin=lciepp, ymax=uciepp), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Ephippia production") +
  theme_bw()

### trait correlation plot
  setkey(epp.ag, clone, gr)
  setkey(male.ag, clone, gr)
  f1 <- merge(epp.ag, male.ag)

  cor.plot <- ggplot(data=f1, aes(x=propmale, y=fillrate, color=gr)) + geom_point() + theme_bw()

### Pool seq results.
    setnames(peaks, "CHROM", "chr")
    setnames(gprime, "CHROM", "chr")
    #setnames(peaks.groups, "CHROM", "chr")
    #setnames(gprime.groups, "CHROM", "chr")

    PA42$hybridchr <- paste("Chr", PA42$PA42chr, "\n", paste("Scaff", "\n", tstrsplit(PA42$chr, "_")[[2]], sep=""), sep="")
    PA42sub <- PA42[, c("chr", "PA42chr", "hybridchr")]
    setkey(PA42sub, chr)
    setkey(peaks, chr)
    setkey(gprime, chr)
    mpeaks <- merge(peaks, PA42sub)
    mgprime <- merge(gprime, PA42sub)

    #peaks <- rbind(peaks, peaks.groups, fill=T)
    #gprime <- rbind(gprime, gprime.groups, fill=T)

    #peaks[,rep:=gsub("NA", "", paste(rep, group, sep=""))]
    #gprime[,rep:=gsub("NA", "", paste(rep, group, sep=""))]

    mgprime[,Gprime.y:=Gprime*ifelse(rep==1, 1, -1)]
    mpeaks[,maxGprime.y:=maxGprime*ifelse(rep==1, 1, -1)]

    mgprime[,pos:=POS]
    setkey(mgprime, hybridchr, pos)
    merge(mgprime,
          data.table(hybridchr=mpeaks$hybridchr,
                     pos=mpeaks$posMaxGprime, key="hybridchr,pos"))[,]

    mgprime$hybridchr <- factor(mgprime$hybridchr, levels=c("Chr1\nScaff\n1931", "Chr2\nScaff\n9198",
      "Chr3\nScaff\n9199", "Chr4\nScaff\n9197", "Chr5\nScaff\n9200", "Chr6\nScaff\n2373",
      "Chr7\nScaff\n7757", "Chr8\nScaff\n6786", "Chr9\nScaff\n1863", "Chr10\nScaff\n2217",
      "Chr11\nScaff\n9201", "Chr12\nScaff\n2158"))

    mpeaks$hybridchr <- factor(mpeaks$hybridchr, levels=c("Chr1\nScaff\n1931", "Chr2\nScaff\n9198",
      "Chr3\nScaff\n9199", "Chr4\nScaff\n9197", "Chr5\nScaff\n9200", "Chr6\nScaff\n2373",
      "Chr7\nScaff\n7757", "Chr8\nScaff\n6786", "Chr9\nScaff\n1863", "Chr10\nScaff\n2217",
      "Chr11\nScaff\n9201", "Chr12\nScaff\n2158"))


    setkey(mgprime, hybridchr)
    setkey(mpeaks, rep, hybridchr)
    poolseq <- ggplot() +
    #geom_vline(data=mpeaks, aes(xintercept=posMaxGprime), color="black") +
    geom_line(data=mgprime, aes(x=pos, y=Gprime.y, group=rep, color=as.factor(rep)), size=1.5) +
    geom_segment(data=mpeaks, aes(x=posMaxGprime, xend=posMaxGprime,
                              y=maxGprime.y+0.15*ifelse(rep==1, 1, -1), yend=5*ifelse(rep==1, 1, -1))) +
    geom_hline(data=mgprime[,list(minG=min(Gprime[qvalue<=0.05])), list(rep)],
           aes(yintercept=minG*ifelse(rep==1, 1, -1))) +
    geom_text(data=mpeaks, aes(x=posMaxGprime, y=5.5*ifelse(rep==1, 1, -1), label=c(1:14))) +
    facet_grid(~hybridchr, scales="free_x", space = "free_x") +
    theme(legend.position = "none") +
    theme_bw() +
    scale_x_continuous(breaks=seq(from=1, to=13, by=2)*10^6, labels=seq(from=1, to=13, by=2)) +
    theme(strip.text.x = element_text(size = 8)) + ylab("Gprime") + xlab("Position (Mb)") +
    labs(color="Replicate")

    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +

### F1 mapping results
  setkey(o.ag.plot, chr)
  setkey(PA42, chr)
  mo.ag.plot <- merge(o.ag.plot, PA42)

  o.ag.perm <- mo.ag.plot[,list(p.aov.thr=mean(p.aov.thr)), list(term, chr, hybridchr, q)]

  o.ag.perm$hybridchr <- factor(o.ag.perm$hybridchr, levels=c("Chr1\nScaff\n1931", "Chr2\nScaff\n9198",
    "Chr3\nScaff\n9199", "Chr4\nScaff\n9197", "Chr5\nScaff\n9200", "Chr6\nScaff\n2373",
    "Chr7\nScaff\n7757", "Chr8\nScaff\n6786", "Chr9\nScaff\n1863", "Chr10\nScaff\n2217",
    "Chr11\nScaff\n9201", "Chr12\nScaff\n2158"))

  mo.ag.plot$hybridchr <- factor(mo.ag.plot$hybridchr, levels=c("Chr1\nScaff\n1931", "Chr2\nScaff\n9198",
    "Chr3\nScaff\n9199", "Chr4\nScaff\n9197", "Chr5\nScaff\n9200", "Chr6\nScaff\n2373",
    "Chr7\nScaff\n7757", "Chr8\nScaff\n6786", "Chr9\nScaff\n1863", "Chr10\nScaff\n2217",
    "Chr11\nScaff\n9201", "Chr12\nScaff\n2158"))

  f1.plot <- ggplot() +
  geom_vline(data=mpeaks, aes(xintercept=posMaxGprime, linetype=as.factor(rep))) +
  geom_line(data=mo.ag.plot, aes(x=pos, y=-log10(p.aov), color=hybridchr), size=1.5) +
  geom_point(data=mo.ag.plot[p.aov<=p.aov.thr], aes(x=pos, y=-log10(p.aov)), size=.5) +
  geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.aov.thr), linetype=as.factor(q))) +
  facet_grid(term~hybridchr, scales="free_x", space = "free_x") +
  theme_bw() + theme(strip.text.x = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(from=1, to=13, by=2)*10^6, labels=seq(from=1, to=13, by=2)) +
  theme(legend.position = "none") + xlab("Position (Mb)")

### overlap plot
  overlap.plot <- ggplot() +
  geom_histogram(data=overlap.perm[perm!=0], aes(z)) +
  geom_vline(data=overlap.perm[perm==0], aes(xintercept=z)) +
  facet_grid(cross~pheno)


### mega plot

phenoplot <- male.plot + epp.plot + emb.plot +
  plot_layout(nrow=1, byrow=TRUE, guides='collect') + plot_annotation(tag_levels = 'A')

mapping <- poolseq + f1.plot +
  plot_layout(nrow-2, byrow=TRUE)

totalCumIndmomwEppMesoOnly <- week7ACtotalind + week7ACmomPEnormWeek + week7ACMomwEppnormWeek + week7ACPropMalenormWeek +
    week7ACEmbWeek + week7ACmomPEnormCumInd + week7ACMomwEppnormCumInd + week7ACPropMalenormCumInd +
    plot_layout(nrow = 2, byrow = TRUE, guides = 'collect') + plot_annotation(tag_levels = 'A')

totalCumIndmomwEppOther <-  Meso + OneLiter + Grace250 + Methyl +
    plot_layout(nrow = 1, byrow = TRUE, guides = 'collect') +
    plot_annotation(tag_levels = list(c("I", "J", "K", "L")))

totaltotal <- totalCumIndmomwEppMesoOnly / totalCumIndmomwEppOther + plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = 'A')


layout <- "
ABC
DDD
EEE
"

bigplot <- male.plot + epp.plot + emb.plot + poolseq + f1.plot +
plot_annotation(tag_levels = 'A', theme=theme(plot.tag = element_text(size = 19))) +
plot_layout(design = layout)

ggsave(bigplot, file="~/big_qtl.png", h=8, w=11)














###
