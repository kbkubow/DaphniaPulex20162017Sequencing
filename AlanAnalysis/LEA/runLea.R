library(LEA)

dap.snmf <- snmf("/scratch/aob2x/daphnia_hwe_sims/lea/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.geno",
    K=c(1:20),
    project = "new",
    repetitions=30,
    CPU = 10,
    entropy=T,
    I=20000,
    iterations = 100000)


save(dap.snmf, file="/scratch/aob2x/daphnia_hwe_sims/lea/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit/dap.snmf.Rdata")
