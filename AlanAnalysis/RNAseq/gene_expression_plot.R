### scp aob2x@rivanna.hpc.virginia.edu:~/ase_geno_phase.star.Rdata ~/.

### libraries
  library(data.table)
  library(ggplot2)

### load output of DEseq2
  load(file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/STAR/DESeq_star.output")

### load ASE data
  load("~/ase_geno.star.Rdata")

### DE subplots
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


  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
    pcaplot

### ASE

  load("~/ase_geno_phase.star.Rdata")

  ase.simple <- ase.geno.phase[ref_dosage==1][allele.x!=allele.y][superclone=="C"][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")]



  qtl_ase <-
  ggplot() +
  geom_hline( yintercept=mean(ase.simple$xCount/ase.simple$totalCount, na.rm=T)) +
  geom_boxplot(data=ase.simple[genes%in%paste("Daphnia0078", c(5:9), sep="")],
                aes(x=clone, y=xCount/totalCount, fill=samp, group=interaction(samp, ref_dosage))) +
  facet_grid(~genes, scales="free_x") + theme_bw()



  ### combined plot
    layout <- "
    AABB
    CDEE
    "

    pcaplot + volcano + magicgene_lfc + magicgene_basemean + geneCounts +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'A') +
    ggtitle("Daphnia00787")














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
