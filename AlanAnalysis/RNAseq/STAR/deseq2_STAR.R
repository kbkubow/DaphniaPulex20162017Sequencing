### scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/rnaseq/featureCounts_trim.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf ~/.
### scp aob2x@rivanna.hpc.virginia.edu:~/featureCounts_STAR.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:~/genes_rd_missing.Rdata ~/.

### libraries
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(Rsubread)
  library("RColorBrewer")
  library(pheatmap)
  library(patchwork)

### load feature counts
  load("~/featureCounts_STAR.Rdata")

### load gene read depth info
  load("~/genes_rd_missing.Rdata")

### SAF/GTF object
  saf <- flattenGTF("~/Daphnia.aed.0.6.gtf",
      method = "merge")

### get sample table
  samps <- as.data.frame(fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable"))
  rownames(samps) <- samps$id
  rownames(samps) <- gsub("pe", ".trim", rownames(samps))
  dds <- DESeqDataSetFromMatrix(countData = fc,
                              colData = samps,
                              design = ~ superclone)


  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]

  dds <- DESeq(dds)
  vsd <- vst(dds, blind=FALSE)

  res <- results(dds, format="DataFrame")
  res

  resLFC <- lfcShrink(dds, coef="superclone_C_vs_A", type="apeglm")

### "best QTL gene"
  plotCounts(dds, gene=which(rownames(dds)=="Daphnia00796"), intgroup="superclone")

  gene_counts <- plotCounts(dds, gene=which(rownames(dds)=="Daphnia00787"), returnData=T, intgroup="superclone", transform=T, normalized=T)
  geneCounts <- ggplot(gene_counts, aes(x=superclone, y=count)) + geom_point() + ylab("Normalized expression - logscale")

### volcano, and LFC density plot
  resLFC$GeneId <- rownames(resLFC)
  res.dt <- as.data.table(resLFC)

  res.dt[GeneId=="Daphnia00787"]
  mean(res.dt$pvalue <= res.dt[GeneId=="Daphnia00787"]$pvalue, na.rm=T)
  mean(abs(res.dt$log2FoldChange) >= abs(res.dt[GeneId=="Daphnia00787"]$log2FoldChange), na.rm=T)

  volcano <- ggplot() +
  geom_point(data=res.dt, aes(x=log2FoldChange, y=-log10(pvalue)), color="red") +
  geom_point(data=res.dt[GeneId=="Daphnia00787"], aes(x=log2FoldChange, y=-log10(pvalue)), color="blue")

  magicgene_lfc <- ggplot(data=res.dt, aes(abs(log2FoldChange))) +
  geom_density() +
  geom_vline(xintercept=abs(res.dt[GeneId=="Daphnia00787"]$log2FoldChange), color="red") +
  geom_text(data=data.frame(x=4, y=1, label=mean(abs(res.dt$log2FoldChange) >= abs(res.dt[GeneId=="Daphnia00787"]$log2FoldChange), na.rm=T)),
            aes(x=x, y=y,
                label=round(label, 3)))

  magicgene_basemean <- ggplot(data=res.dt, aes(log10(baseMean))) +
  geom_density() +
  geom_vline(xintercept=log10(res.dt[GeneId=="Daphnia00787"]$baseMean), color="red") +
  geom_text(data=data.frame(x=4, y=1, label=mean(log10(res.dt$baseMean) >= log10(res.dt[GeneId=="Daphnia00787"]$baseMean), na.rm=T)),
            aes(x=x, y=y,
                label=round(label, 3)))

### PCA analysis
  pcaData <- plotPCA(vsd, intgroup=c("superclone", "clone"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
    pcaplot


### volcano combined with genes.ag
  genes.ag.ag <- genes.ag[,list(scA=mean(dp.norm.mu[superclone=="A"], na.rm=T) , scB=mean(dp.norm.mu[superclone=="C"], na.rm=T),
                                  missingA=mean(missing.rate[superclone=="A"], na.rm=T) , missingB=mean(missing.rate[superclone=="C"], na.rm=T)), list(GeneId=GeneID)]
  res.gene <- merge(res.dt, genes.ag.ag, by="GeneId")
  res.gene[,delta:=scA-scB]
  res.gene[,quanLFC:=rank(log2FoldChange)]
  res.gene[quanLFC==1]

  ggplot(data=res.gene, aes(y=delta, x=log2FoldChange)) + geom_point()

  ggplot(data=res.gene, aes(x=scA, y=scB, color=log2FoldChange)) + geom_point()



  fisher.test(table(res.gene$missingB==1 & res.gene$missingA<1, res.gene$log2FoldChange< -4))



  summary(lm(log2FoldChange~scA*scB, res.gene))



### combined plot
  layout <- "
  AABB
  CDEE
  "

  pcaplot + volcano + magicgene_lfc + magicgene_basemean + geneCounts +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') +
  ggtitle("Daphnia00787")


### combine with QTL
  load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")
  peaks.dt <- data.table(qtl=c(1:14), Chr=peaks$CHROM, Start=peaks$posPeakDeltaSNP-50000, End=peaks$posPeakDeltaSNP+50000)

  saf.dt <- as.data.table(saf)
  setkey(saf.dt, GeneID)
  saf.dt<-unique(saf.dt, by = key(saf.dt), fromLast = TRUE)

  setkey(peaks.dt, Chr, Start, End)
  setkey(saf.dt, Chr, Start, End)

  qtl.genes <- foverlaps(peaks.dt, saf.dt)

  setnames(qtl.genes, "GeneID", "GeneId")
  setkey(qtl.genes, "GeneId")

  resLFC$GeneId <- row.names(resLFC)
  resLFC <- as.data.table(resLFC)

  res.qtl <- merge(resLFC, qtl.genes[,c("GeneId", "qtl"),with=F], by="GeneId", all.x=T)
  setkey(res.qtl, "GeneId")

  res.qtl[!is.na(qtl)][order(padj)]
  setnames(saf.dt, "GeneID", "GeneId")
  res.qtl <- merge(res.qtl, saf.dt, by="GeneId")

  res.qtl.ag <- res.qtl[,list(.N), Chr]

  dim(res.qtl)
  res.qtl <- res.qtl[Chr%in%res.qtl.ag[1:12]$Chr]
  dim(res.qtl)

  res.qtl[,mid:=Start/2+End/2]
  res.qtl[,quan:=rank(padj)/(length(padj)-1)]

  res.qtl.ag <- res.qtl[,list(frac=mean(quan<=.1),
                              TT=sum(quan<=.1),
                              expected=.1*.N,
                              .N,
                              p=binom.test(sum(quan<=.1), sum(quan>.1), .1)$p.value), list(qtl)]

  res.qtl[qtl==10][quan<=.1]

### plot
  res.qtl[,mid:=Start/2+End/2]
  res.qtl[,quan:=rank(padj)/(length(padj)-1)]

  ggplot(data=res.qtl[order(quan, na.last=F)], aes(x=abs(log2FoldChange), y=pvalue, color=!is.na(qtl),
        shape=padj<1e-20)) +
  geom_point(size=1, alpha=.75)


  ggplot(data=res.qtl[order(qtl, na.last=F)], aes(x=log2FoldChange, y=-log10(pvalue), color=!is.na(qtl))) + geom_point()


    ggplot(data=res.qtl[order(qtl, na.last=F)], aes(x=mid, y=-log10(padj), color=!is.na(qtl),
          shape=padj<1e-10)) +
    geom_point(size=.75, alpha=.75) + facet_wrap(~Chr, nrow=1, scales="free_x")



  scp aob2x@rivanna.hpc.virginia.edu:~/res_deseq.Rdata ~/.


  library(data.table)
  library(ggplot2)
  load("~/res_deseq.Rdata")

  ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj))) + geom_point()
