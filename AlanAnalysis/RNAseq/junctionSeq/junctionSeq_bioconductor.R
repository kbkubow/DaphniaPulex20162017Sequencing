# module load gcc/7.1.0
# ~/R-3.2.2/bin/R

#LoadR JunctionSeq:
  library("JunctionSeq")

#The sample decoder:
  decoder <- read.table("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/decoder.txt", header=T,stringsAsFactors=F)

#The count files:
  countFiles <- paste0("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/",
                      decoder$sample.ID,
                      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

#Run the analysis:
  jscs <- runJunctionSeqAnalyses(
    sample.files = countFiles,
    sample.names = decoder$sample.ID,
    condition = decoder$group.ID,
    use.novel.junctions=T,
    use.multigene.aggregates=T,
    flat.gff.file = "/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/withNovel.forJunctionSeq.gff.gz",
    verbose=TRUE, debug.mode = TRUE, nCores=20, method.GLM="simpleML", method.dispFit="mean")

# save
  save(jscs, file="~/jscs.Rdata")

# plot
  buildAllPlots(jscs=jscs,
    outfile.prefix = "~/jcsc_plots/",
    use.plotting.device = "png", FDR.threshold = 0.01, writeHTMLresults=T
  )
