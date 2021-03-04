#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R



### library
  library(data.table)
  library(foreach)
  library(stringr)

### library size data

  samps <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable")
  samps[,unique.ID:=gsub("pe.bam", "", id)]

  nreads <- foreach(samp.i=samps$SampleName, .combine="rbind")%do%{
    #samp.i <- samps[1]$unique.ID
    nReads <- as.numeric(str_match(
              system(paste("grep  \"Number of input reads\" ",
                "/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/", samp.i, "_star_testLog.final.out", sep=""), intern=T),
                "[0-9]{1,}"))
    data.table(samp=samp.i, nReads=nReads)
  }

### load data
  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out", full.names=T)

  genebody <- foreach(fn.i=fn, .errorhandling="remove")%do%{
    #fn.i <- fn[1]
    tmp <- fread(paste(fn.i, "QC.geneBodyCoverage.genewise.txt.gz", sep="/"), header=T)

    gbl <- melt(tmp, "GENE_ID")
    gbl.ag <- gbl[,list(aveExp=mean(value)), list(GENE_ID)]
    gbl <- merge(gbl, gbl.ag, by="GENE_ID")
    gbl[,samp:=last(tstrsplit(fn.i, "/"))]
    gbl <- merge(gbl, nreads, by="samp")

  }
  genebody <- rbindlist(genebody)
  genebody <- merge(genebody, samps, by.x="samp", by.y="SampleName")

  save(genebody, file="~/genebody.Rdata")

# scp aob2x@rivanna.hpc.virginia.edu:~/genebody.Rdata ~/.
  library(data.table)
  library(ggplot2)

  load("~/genebody.Rdata")

  genebody[,exp.norm:=(value/nReads)/(aveExp/nReads)]
  genebody[,var:=as.numeric(as.character(variable))]

  genebody.ag <- genebody[,list(mu=mean(exp.norm, na.rm=T), sd=sd(exp.norm, na.rm=T)), list(var=var)]

  ggplot() +
  geom_ribbon(data=genebody.ag, aes(x=var, ymax=mu-2*sd, ymin=mu+2*sd), alpha=.3) +
  geom_line(data=genebody.ag, aes(x=var, y=mu)) +
  geom_line(data=genebody[GENE_ID=="Daphnia00787"], aes(x=var, y=exp.norm, group=samp, color=superclone)) +
  ylab("Normalized gene expression") +
  xlab("Distance along genebody 5'->3'") +
  theme_bw()
