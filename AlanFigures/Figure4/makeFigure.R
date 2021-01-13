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
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")

### male proportion plot
  male.ag <- male[!is.na(gr),list(propmale=sum(Males)/sum(NewTotal),
                                  N=sum(NewTotal)),
                    list(clone, gr)]
  male.ag[,se:=sqrt((propmale*(1-propmale))/N)]

  male.ag[,lci:=propmale-1.96*se]
  male.ag[,uci:=propmale+1.96*se]

  male.plot <- ggplot(data=male.ag, aes(x=gr, y=propmale, color=gr)) +
  geom_linerange(aes(ymin=lci, ymax=uci), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Male proportion") +theme_bw()

### Fill Rate
  epp.ag <- epp[!is.na(gr),list(fillrate=sum(fill*TotalEppB, na.rm=T)/sum(TotalEppB),
                                N=sum(TotalEppB)),
                  list(clone, gr)]

  epp.ag[,se:=sqrt((fillrate*(1-fillrate))/N)]
  epp.ag[,lci:=fillrate-1.96*se]
  epp.ag[,uci:=fillrate+1.96*se]

  epp.plot <- ggplot(data=epp.ag, aes(x=gr, y=fillrate, color=gr)) +
  geom_linerange(aes(ymin=lci, ymax=uci), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Ephippial fill rate") +
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

  ### we are takign the average just because of a consequence of how the data are output. e.g., `o[id==111][term=="male"][perm==0]`
  o.ag <- o[, list(p.z=mean(p.z, na.rm=T), chisq=mean(chisq), p.aov=mean(p.aov, na.rm=T)), list(term, perm, id, chr, pos)]
  o.ag.ag <- o.ag[,list(pr=1-mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id, chr, pos)]
  o.ag.ag[pr==0, pr:=1/201]

  ### ANOVA output
  o.ag.perm <- o.ag[perm!=0, list(p.aov.thr=c(quantile(p.aov, .01), quantile(p.aov, .05)),
                                  q=c(.01, .05)), list(term, chr)]

  setkey(o.ag.perm, term, chr)
  setkey(o.ag, term, chr)
  o.ag.plot <- merge(o.ag[perm==0], o.ag.perm[q==0.01], allow.cartesian=T)


  f1.plot <- ggplot() +
  geom_vline(data=peaks, aes(xintercept=posMaxGprime, linetype=as.factor(rep))) +
  geom_line(data=o.ag.plot, aes(x=pos, y=-log10(p.aov), color=chr), size=1.5) +
  geom_point(data=o.ag.plot[p.aov<=p.aov.thr], aes(x=pos, y=-log10(p.aov)), size=.5) +
  geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.aov.thr), linetype=as.factor(q))) +
  facet_grid(term~chr, scales="free_x") +
  theme(legend.position = "none")


### mega plot
  (male.plot + epp.plot + cor.plot) / poolseq / f1.plot
















###
