library(data.table)
library(foreach)

fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/coverage/", "coveragePos.delim")


scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/rnaseq/coverage/d8_349_2_Daphnia00787.coveragePos.delim ~/.

library(data.table)
library(ggplot2)



dat <- fread("~/d8_349_2_Daphnia00787.coveragePos.delim")

ggplot(data=dat, aes(x=V3, y=V4)) +
geom_line() +
geom_hline(yintercept=0)


plot(V4~V3, dat)
abline(h=0)
