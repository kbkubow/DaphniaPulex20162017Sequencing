scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/rnaseq/quorts.tar.gz ~/.


install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz",
                   repos=NULL,
                   type="source");
library(QoRTs)
library(data.table)

dec <- fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable")
dec[,unique.ID:=gsub("pe.bam", "", id)]

decoder <- data.table(unique.ID=gsub("pe.bam", "", dec$id),
                      lane.ID="UNKNOWN",
                      group.ID=dec$superclone,
                      sampleID=dec$rep)
decoder[,qc.data.dir:=paste("/Users/alanbergland/qorts_out", unique.ID, sep="/")]


res <- read.qc.results.data(infile.dir="/",
  decoder=as.data.frame(decoder)[1,],
  calc.DESeq2 = FALSE,
  calc.edgeR = FALSE,
  autodetectMissingSamples = TRUE, skip.files = c())
plotter.basic <- build.plotter.basic(res)
makePlot.genebody.coverage(plotter.basic)

makePlot.genebody(plotter.basic,
                  geneset = c("Overall","90-100","75-90","50-75","0-50"),
                  avgMethod = c("TotalCounts", "AvgPercentile"),
                  plot.medians,
                  plot.means = TRUE,
                  debugMode,
                  singleEndMode,
                  ... )
