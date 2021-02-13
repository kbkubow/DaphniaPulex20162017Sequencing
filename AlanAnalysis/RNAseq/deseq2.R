### scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/rnaseq/featureCounts_trim.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf ~/.

### libraries
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(Rsubread)
  library("RColorBrewer")
  library(pheatmap)

### load feature counts
  load("~/featureCounts_trim.Rdata")

### SAF/GTF object
  saf <- flattenGTF("~/Daphnia.aed.0.6.gtf",
      method = "merge")

### get sample table
  samps <- as.data.frame(fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable"))
  rownames(samps) <- samps$id
  rownames(samps) <- gsub("pe", ".trim", rownames(samps))
  dds <- DESeqDataSetFromMatrix(countData = fc$counts,
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

### "best gene"
  plotCounts(dds, gene=which.min(res$padj), intgroup="superclone")
  dds <- estimateSizeFactors(dds)

### PCA analysis
  pcaData <- plotPCA(vsd, intgroup=c("superclone", "clone"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()



### clustering
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)

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

  res$GeneId <- row.names(res)
  res <- as.data.table(res)

  res.qtl <- merge(res, qtl.genes[,c("GeneId", "qtl"),with=F], by="GeneId", all.x=T)
  setkey(res, "GeneId")

  res[!is.na(qtl)][order(padj)]

  scp aob2x@rivanna.hpc.virginia.edu:~/res_deseq.Rdata ~/.


  library(data.table)
  library(ggplot2)
  load("~/res_deseq.Rdata")

  ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj))) + geom_point()
