### scp aob2x@rivanna.hpc.virginia.edu:~/ase_geno_phase.star.Rdata ~/.

### libraries
  library(data.table)
  library(ggplot2)
  library(DESeq2)
  library(tidyverse)
  library(foreach)

### load output of DEseq2
  load(file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/STAR/DESeq_star.output")

### load ASE data
  load("~/ase_geno.star.Rdata")

  targetGenes <- paste("Daphnia0078", c(6:9), sep="")


### load nominal read depth

  load("~/coverage_pos_rna.Rdata")
  setnames(dat, names(dat), c("samp", "chrLen", "readsMapped", "chr", "pos", "rd"))
  pcaData <- as.data.table(pcaData)
  pcaData[,samp:=gsub(".trim.bam", "", name)]
  dat <- merge(dat, pcaData, by="samp")

  dat[,rd.norm:=rd/readsMapped]
  dat[pos>=5191562-500 & pos<=5191562, gene:="Daphnia00787_3"]
  dat[pos>=5191562 & pos<=5204101, gene:="Daphnia00787"]
  dat[pos>=5204101 & pos<=5204101+5000, gene:="Daphnia00787_5"]
  dat[is.na(gene), gene:="other"]

  norm.rd <- ggplot(data=dat[!is.na(gene)], aes(x=pos, y=rd.norm, group=interaction(samp, gene), color=superclone, linetype=gene)) +
  geom_line() +
  geom_hline(yintercept=0)

  real.rd <- ggplot(data=dat[!is.na(gene)], aes(x=pos, y=rd, group=interaction(samp, gene), color=superclone, linetype=gene)) +
  geom_line() +
  geom_hline(yintercept=0)

  norm.rd + real.rd
### DE subplots
  volcano <- ggplot() +
  geom_point(data=res.dt, aes(x=log2FoldChange, y=-log10(pvalue)), color="red") +
  geom_point(data=res.dt[GeneId%in%targetGenes], aes(x=log2FoldChange, y=-log10(pvalue)), color="blue")

  magicgene_lfc <- ggplot(data=res.dt, aes(abs(log2FoldChange))) +
  geom_density() +
  geom_vline(xintercept=abs(res.dt[GeneId%in%targetGenes]$log2FoldChange), color="red") +
  geom_text(data=data.frame(x=6, y=seq(from=2, to=1, length.out=4),
                                      label=sapply(targetGenes, function(x) mean(abs(res.dt$log2FoldChange) >= abs(res.dt[GeneId==x]$log2FoldChange), na.rm=T))%>%round(., 3)%>%paste(targetGenes, ., sep=": ")),
            aes(x=x, y=y,
                label=label))

  magicgene_pvalue <- ggplot(data=res.dt, aes(log10(-log10(pvalue)))) +
  geom_density() +
  geom_vline(xintercept=abs(res.dt[GeneId%in%targetGenes]$log2FoldChange), color="red") +
  geom_text(data=data.frame(x=-2.5, y=seq(from=2, to=1, length.out=4),
                                      label=sapply(targetGenes, function(x) mean((res.dt$pvalue) <= (res.dt[GeneId==x]$pvalue), na.rm=T))%>%round(., 3)%>%paste(targetGenes, ., sep=": ")),
            aes(x=x, y=y,
                label=label))


  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
    pcaplot


  gene_counts <- foreach(i=6:9, .combine="rbind")%do%{
    gene_counts <- as.data.table(plotCounts(dds, gene=which(rownames(dds)==paste("Daphnia0078", i, sep="")), returnData=T, intgroup="superclone", transform=T, normalized=T))
    gene_counts[,gene:=paste("Daphnia0078", i, sep="")]
  }

  geneCounts <- ggplot(gene_counts, aes(x=superclone, y=count)) + geom_point() + ylab("Normalized expression - logscale") + facet_wrap(~gene, scales="free_y", nrow=1)



### ASE

  load("~/ase_geno_phase.star.Rdata")

  ase.simple <- ase.geno.phase[ref_dosage==1][allele.x!=allele.y][superclone=="C"][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")]



  qtl_ase <-
  ggplot() +
  geom_hline( yintercept=mean(ase.simple$xCount/ase.simple$totalCount, na.rm=T)) +
  geom_boxplot(data=ase.simple[genes%in%targetGenes],
                aes(x=clone, y=xCount/totalCount, fill=samp, group=interaction(samp, ref_dosage))) +
  facet_grid(~genes, scales="free_x") + theme_bw()


### read depth

  dep.plot <- ggplot(data=ase.geno.phase[genes%in%c("Daphnia00787")], aes(x=pos, y=totalCount, group=samp, color=superclone)) +
    geom_line() + theme(legend.position="bottom")

  ase.plot <- ggplot(data=ase.geno.phase[genes%in%c("Daphnia00787")][ref_dosage==1][superclone=="C"], aes(x=pos, y=xCount/totalCount, group=samp, color=class)) +
    geom_point() + theme(legend.position="bottom")

  biplot.plot <- ggplot(data=ase.geno.phase[genes%in%c("Daphnia00787")][ref_dosage==1][superclone=="C"], aes(x=totalCount, y=xCount/totalCount, group=samp, color=superclone)) +
    geom_point() + theme(legend.position="bottom")

### pnps distributions
  annotation.counts <- ase.geno.phase[,list(.N), list(superclone, samp, class, genes)]
  annotation.counts.ag <- annotation.counts[,list(pn=N[class=="missense_variant"], ps=N[class=="synonymous_variant"]), list(superclone, samp, genes)]

  ggplot() +
  geom_histogram(data=annotation.counts.ag[pn!=0 & ps!=0 & !is.na(pn) & !is.na(ps)], aes(log2(pn/ps))) +
  geom_vline(data=annotation.counts.ag[genes=="Daphnia00787"], aes(xintercept=log2(pn/ps)), color="red")

  ### combined plot
    layout <- "
    AAABBB
    CCCEEE
    FFFFFF
    GGGGGG
    HHHIIJ
    HHHIIJ
    "

    big.plot <-
    pcaplot + volcano +
    magicgene_lfc + magicgene_pvalue +
    geneCounts +
    qtl_ase +
    dep.plot + ase.plot + biplot.plot +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A')

    ggsave(big.plot, file="~/expression_plot.pdf")










  ### ASE subplots
    ase.geno.ag <- ase.geno[class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")][,
                            list(ai=mean(refCount/totalCount, na.rm=T), .N), list(gene=genes, geno=samp, ref_dosage)]

    ase.geno[ref_dosage==1,p1:=pbinom(refCount, totalCount, .5)]
    ase.geno[ref_dosage==1,p2:=1-pbinom(refCount, totalCount, .5)]

    ggplot(ase.geno[genes=="Daphnia00787"][ref_dosage==1][p1<.005 | p2<.005], aes(x=pos, y=refCount/totalCount, color=superclone)) + geom_point()

    ggplot(data=ase.geno.ag[ref_dosage==1],
          aes(ai, group=geno, color=geno)) +
    geom_density()

  ###
