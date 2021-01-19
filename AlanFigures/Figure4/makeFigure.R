### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(patchwork)

### setwd
  setwd("/Users/alanbergland/Documents/GitHub")

### load F1 cross data (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")

### load in PoolSeq data (made by `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLSeqR.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")

### load in F1 mapping data
  #load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")
  #load("~/lme4qtl_output.AxC.long.Rdata")
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.SUMMARIZED.Rdata")

### load in overlap tests
  load("~/overlap_perm.Rdata")

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

    #peaks <- rbind(peaks, peaks.groups, fill=T)
    #gprime <- rbind(gprime, gprime.groups, fill=T)

    #peaks[,rep:=gsub("NA", "", paste(rep, group, sep=""))]
    #gprime[,rep:=gsub("NA", "", paste(rep, group, sep=""))]

    gprime[,Gprime.y:=Gprime*ifelse(rep==1, 1, -1)]
    peaks[,maxGprime.y:=maxGprime*ifelse(rep==1, 1, -1)]

    gprime[,pos:=POS]
    setkey(gprime, chr, pos)
    merge(gprime,
          data.table(chr=peaks$chr,
                     pos=peaks$posMaxGprime, key="chr,pos"))[,]


    setkey(gprime, chr)
    poolseq <- ggplot() +
    #geom_vline(data=peaks, aes(xintercept=posMaxGprime), color="black") +
    geom_line(data=gprime, aes(x=pos, y=Gprime.y, group=rep, color=as.factor(rep)), size=1.5) +
    geom_segment(data=peaks, aes(x=posMaxGprime, xend=posMaxGprime,
                                  y=maxGprime.y+0.15*ifelse(rep==1, 1, -1), yend=5*ifelse(rep==1, 1, -1))) +
    geom_hline(data=gprime[,list(minG=min(Gprime[qvalue<=0.05])), list(rep)],
               aes(yintercept=minG*ifelse(rep==1, 1, -1))) +
    geom_text(data=peaks, aes(x=posMaxGprime, y=5.5*ifelse(rep==1, 1, -1), label=c(1:14))) +
    facet_grid(~paste("Scaffold_", tstrsplit(chr, "_")[[2]], sep=""), scales="free_x", space = "free_x") +
    theme(legend.position = "none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    scale_x_continuous(breaks=seq(from=1, to=9, by=2)*10^6, labels=seq(from=1, to=9, by=2))

### F1 mapping results
  o.ag.perm <- o.ag.plot[,list(p.aov.thr=mean(p.aov.thr)), list(term, chr, q)]

  f1.plot <- ggplot() +
  geom_vline(data=peaks, aes(xintercept=posMaxGprime, linetype=as.factor(rep))) +
  geom_line(data=o.ag.plot, aes(x=pos, y=-log10(p.aov), color=chr), size=1.5) +
  geom_point(data=o.ag.plot[p.aov<=p.aov.thr], aes(x=pos, y=-log10(p.aov)), size=.5) +
  geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.aov.thr), linetype=as.factor(q))) +
  facet_grid(term~chr, scales="free_x") +
  theme(legend.position = "none")

### overlap plot
  overlap.plot <- ggplot() +
  geom_histogram(data=overlap.perm[perm!=0], aes(z)) +
  geom_vline(data=overlap.perm[perm==0], aes(xintercept=z)) +
  facet_grid(cross~pheno)


### mega plot
layout <- "
ABCD
EEEE
FFFF
"

bigplot <- male.plot + epp.plot + cor.plot + overlap.plot + poolseq + f1.plot +
plot_annotation(tag_levels = 'A', theme=theme(plot.tag = element_text(size = 19))) +
plot_layout(design = layout)

ggsave(bigplot, file="~/big_qtl.png", h=8, w=11)














###
