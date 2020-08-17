library(LEA)
dap.snmf <- snmf("/scratch/aob2x/daphnia_hwe_sims/daps4lea.geno",
    K=c(1:10),
    project = "continue",
    repetitions = 1, CPU = 1,
    alpha = 10, tolerance = 0.00001, entropy = FALSE, percentage = 0.05, I=10000, iterations = 200, ploidy = 2, seed = -1)


save(dap.snmf, file="/scratch/aob2x/daphnia_hwe_sims/dap.snmf.Rdata")
