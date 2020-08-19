library(LEA)

dap.snmf <- snmf("/scratch/aob2x/daphnia_hwe_sims/daps4lea.geno",
    K=c(1:10),
    project = "continue",
    repetitions=10,
    CPU = 10,
    entropy=T,
    I=20000)


save(dap.snmf, file="/scratch/aob2x/daphnia_hwe_sims/dap.snmf.Rdata")
