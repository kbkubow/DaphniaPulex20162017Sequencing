scp aob2x@rivanna.hpc.virginia.edu:/home/aob2x/cdlo_250K.boot.manhattan.Rdata ~/.

library(ggplot2)
library(data.table)

load("~/cdlo_250K.boot.manhattan.Rdata")
load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")
setnames(peaks, "CHROM", "chr")


temp <- cdl.o.manhattanPlot.ag[!is.na(pond.group)][,list(rat=cd_mean[pond.group=="A-C"]/cd_mean[pond.group=="DWT-DWT"]), list(chr, mid)]

mp_small <- ggplot(data=cdl.o.manhattanPlot.ag, aes(x=mid, y=log2(rat), color=chr)) +
geom_vline(data=peaks, aes(xintercept=posMaxGprime)) +
geom_line(size=1) +
facet_grid(~chr, scales="free")

ggsave(mp_small, file="~/mp_small.boot.png", height=10, w=20)


mp_small <- ggplot(data=cdl.o.manhattanPlot.ag[!is.na(pond.group)], aes(x=mid, y=cd_max, color=chr)) +
geom_vline(data=peaks, aes(xintercept=posMaxGprime)) +
geom_line(size=1) +
facet_grid(sp.group+pond.group~chr, scales="free")

ggsave(mp_small, file="~/mp_small.boot.png", height=10, w=20)
