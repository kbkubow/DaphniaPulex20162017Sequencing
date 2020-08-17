library(LEA)
dap.snmf <- snmf("/scratch/aob2x/daphnia_hwe_sims/daps4lea.geno",
    K=c(1:10),
    project = "new",
    CPU = 10, entropy=T)


save(dap.snmf, file="/scratch/aob2x/daphnia_hwe_sims/dap.snmf.Rdata")
