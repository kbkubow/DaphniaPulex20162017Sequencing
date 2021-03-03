/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.flatgff



suppressPackageStartupMessages(library(DEXSeq));
decoder <- read.table(
"inputData/annoFiles/decoder.bySample.txt", header=T,stringsAsFactors=F);
23
directory <- "outputData/countTables/";
countFiles <- paste0( directory,
decoder$sample.ID,
"/QC.exonCounts.formatted.for.DEXSeq.txt.gz"
);
dexseq.gff <- "outputData/forDEXSeq.gff.gz";
sampleTable <- data.frame( row.names = decoder$sample.ID, condition = decoder$group.ID
);
dxd <- DEXSeqDataSetFromHTSeq( countFiles,
);
