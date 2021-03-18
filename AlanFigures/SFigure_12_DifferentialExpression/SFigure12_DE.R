### libraries
  library(data.table)
  library(ggplot2)
  library(DESeq2)
  library(tidyverse)
  library(gdata)

### load general DE objects
  setwd("/Users/alanbergland/Documents/GitHub/")

  load(file="DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/differential_expression_supplement.Rdata")
  load(file="DaphniaPulex20162017Sequencing/AlanFigures/Figure4/differentialExpression_withNames.Rdata")

  #### some data stuff
    ### PCA data
      pcaData <- as.data.table(pcaData)
      pcaData[,clone:=gsub("d8_", "", clone)]
      pcaData[,rep:=tstrsplit(name, "_")[[3]]%>%gsub(".trim.bam", "", .)]

### load in annotation sheet
  go <- read.xls("DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/Daphnia_annotation_PANTHER.xls")
  go <- as.data.table(go)


### go analysis
  go.dt <- go[,c("qseqid", "GO"), with=F]
  go.dt <- go.dt[GO!=""]
  go.dt[,gene:=tstrsplit(qseqid, "-")[[1]]]
  go.dt.ag <- go.dt[,list(GO=paste(GO, sep=";")), list(gene)]

  geneGOid <- foreach(i=1:dim(go.dt.ag)[1])%do%{
    strsplit(go.dt.ag[i]$GO, ";")[[1]]
  }
  names(geneGOid) <- go.dt.ag$gene

  allGenes <- dec[!is.na(qval)]$qval
  names(allGenes) <- dec[!is.na(qval)]$GeneID

  geneNames <- names(geneGOid)
  head(geneNames)

  myInterestingGenes <- dec[!is.na(qval)][qval<.05]$GeneID
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  str(geneList)


  GOdata <- new("topGOdata",
                ontology="MF",
                allGenes=geneList,
                annot=annFUN.gene2GO,
                gene2GO=geneGOid)
  graph(GOdata)
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")


  allRes <- GenTable(GOdata, classic = resultFis,
                  orderBy = "weight", ranksOf = "classic", topNodes = 20)


### plot components
  ### PCA
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      theme_bw()

  ###copy number DE plot
    cn.plot <- ggplot(data=dec[(cnA-cnB)%in%c(-2, -1, 0, 1, 2)][cnA%in%c(0,1,2) & cnB%in%c(0,1,2)], aes(x=as.factor(I(cnB-cnA)), y=log2FoldChange)) +
    geom_boxplot() + theme_bw() +
    ylab("log2(Fold Change): C / A") +
    xlab("(Copy Number C) - (Copy Number A)")




















Daphnia08075
plotCounts(dds, gene=which(rownames(dds)=="Daphnia08075"), returnData=F, intgroup="superclone", transform=T, normalized=T)

### gene expression plots

  gene_counts <- foreach(i=dec[qLFC.noCN.goodChr>.9 & !is.na(old_QTL_ID)]$GeneID, .combine="rbind")%do%{
    gene_counts <- as.data.table(plotCounts(dds, gene=which(rownames(dds)==i), returnData=T, intgroup="superclone", transform=T, normalized=T))
    gene_counts[,gene:=i]
  }
  gene_counts <- merge(gene_counts, dec, by.x="gene", by.y="GeneID")

  table(gene_counts$final_QTL_ID)


  ### QTL 5: 3 genes
  qtl_5_de <-
    ggplot(gene_counts[final_QTL_ID==5],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_grid(PA42qtl~gene, scales="free") + theme_bw()

  ### QTL 10: 3 genes
  qtl_8_de <-
    ggplot(gene_counts[final_QTL_ID==10],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_grid(PA42qtl~gene, scales="free") + theme_bw()

  ### QTL 12: 4 genes
  qtl_12_de <-
    ggplot(gene_counts[final_QTL_ID==12],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_wrap(PA42qtl~gene, scales="free", nrow=1) + theme_bw()

  ### QTL 7: 1 genes
  qtl_7_de <-
    ggplot(gene_counts[final_QTL_ID==7],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_grid(PA42qtl~gene, scales="free") + theme_bw()

  qtl_9_de <-
    ggplot(gene_counts[final_QTL_ID==9],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_grid(PA42qtl~gene, scales="free") + theme_bw()


qtl_5_de + plot_spacer() + plot_spacer() + plot_spacer() / qtl_8_de / qtl_9_de / qtl_7_de / qtl_9_de


### plot components
  ### PCA
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      theme_bw()

  ###
