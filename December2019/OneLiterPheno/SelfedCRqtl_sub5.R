ijob -A berglandlab -p standard -t 05:00:00 -N 1 -n 1 --mem 200000
module load  gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3

#!/usr/bin/env Rscript

library(qtl)
setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/selfC")
selfC <- read.cross("csv","","selfedCforrqtlwithpheno_sub5.csv")
summary(selfC)

    F2 intercross

    No. individuals:    23

    No. phenotypes:     3
    Percent phenotyped: 100 95.7 100

    No. chromosomes:    11
        Autosomes:      Scaffold_1863_HRSCAF_2081 Scaffold_1931_HRSCAF_2197
                        Scaffold_2158_HRSCAF_2565 Scaffold_2217_HRSCAF_2652
                        Scaffold_2373_HRSCAF_2879 Scaffold_6786_HRSCAF_7541
                        Scaffold_7757_HRSCAF_8726 Scaffold_9197_HRSCAF_10753
                        Scaffold_9199_HRSCAF_10755 Scaffold_9200_HRSCAF_10757
                        Scaffold_9201_HRSCAF_10758

    Total markers:      1450
    No. markers:        30 160 120 460 130 80 120 40 130 100 80
    Percent genotyped:  88.6
    Genotypes (%):      AA:25.0  AB:54.0  BB:20.9  not BB:0.0  not AA:0.0

plotMissing(selfC)
par(mfrow=c(1,2), las=1)
plot(ntyped(selfC), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(selfC, "mar"), ylab="No. typed individuals", main="No. genotypes by marker")

cg <- comparegeno(selfC)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

# No matching genotypes, which is as it should be
#Dropping individuals with segregation distortion.
gt <- geno.table(selfC)
gt[gt$P.value < 0.05/totmar(selfC),]
todrop <- rownames(gt[gt$P.value < 1e-10,])
# No markers below threshold to drop

selfC <- est.rf(selfC)
plotMissing(selfC)
plotMap(selfC)
plotRF(selfC)

plotRF(selfC, chr="Scaffold_1931_HRSCAF_2197")

summaryMap(selfC)

selfC <- switch.order(selfC, chr="Scaffold_1931_HRSCAF_2197", c(1:33, 64:65, 34:63, 66:85), error.prob=0.005)
selfC <- switch.order(selfC, chr="Scaffold_1931_HRSCAF_2197", c(1:32, 34:37, 33, 38:85), error.prob=0.005)

# Let's just leave this for now.

plotRF(selfC, chr="Scaffold_1863_HRSCAF_2081")

# Looks fine

plotRF(selfC, chr="Scaffold_2158_HRSCAF_2565")

selfC <- switch.order(selfC, chr="Scaffold_2158_HRSCAF_2565", c(1:45, 53, 46:52, 54:60), error.prob=0.005)

# OK

plotRF(selfC, chr="Scaffold_2217_HRSCAF_2652")

#Leaving for now

plotRF(selfC, chr="Scaffold_2373_HRSCAF_2879")

#Ok good

plotRF(selfC, chr="Scaffold_6786_HRSCAF_7541")

#Fine

plotRF(selfC, chr="Scaffold_7757_HRSCAF_8726")

selfC <- switch.order(selfC, chr="Scaffold_7757_HRSCAF_8726", c(1:56, 61:70, 57:60), error.prob=0.005)
selfC <- switch.order(selfC, chr="Scaffold_7757_HRSCAF_8726", c(1:40, 51:55, 41:50, 56:70), error.prob=0.005)

#OK-ish

plotRF(selfC, chr="Scaffold_9197_HRSCAF_10753")

#Fine

plotRF(selfC, chr="Scaffold_9199_HRSCAF_10755")

selfC <- switch.order(selfC, chr="Scaffold_9199_HRSCAF_10755", c(1:45, 81:85, 46:80), error.prob=0.005)
selfC <- switch.order(selfC, chr="Scaffold_9199_HRSCAF_10755", c(1:25, 31:45, 26:30, 46:85), error.prob=0.005)
selfC <- switch.order(selfC, chr="Scaffold_9199_HRSCAF_10755", c(1:16, 21:25, 17:20, 26:85), error.prob=0.005)

# OK

plotRF(selfC, chr="Scaffold_9200_HRSCAF_10757")

# Fine

plotRF(selfC, chr="Scaffold_9201_HRSCAF_10758")

selfC <- switch.order(selfC, chr="Scaffold_9201_HRSCAF_10758", c(1:48, 60, 49:60), error.prob=0.005)
selfC <- switch.order(selfC, chr="Scaffold_9201_HRSCAF_10758", c(1:49, 60, 50:59), error.prob=0.005)

# OK

selfC <- est.rf(selfC)
plotRF(selfC)
plotMap(selfC)

selfC <- calc.genoprob(selfC,step=2.5,error.prob=0.01)

save.image(file="selfC20200709_sub_A.RData")



#!/usr/bin/env Rscript

library(qtl)
setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/selfC")

load("selfC20200709_sub_A.RData")
selfC <- sim.geno(selfC, step=2.5, n.draws=64, error.prob=0.01)

save.image(file="selfC20200708_D.RData")
save.image(file="selfC20200708_E.RData")

out.em <-scanone(hyper)
out.imp <-scanone(hyper,method="imp")
outepp.em <- scanone(selfC, pheno.col=1)
outeppimp.em <- scanone(selfC, pheno.col=1, method="imp")
outemb.em <- scanone(selfC, pheno.col=2)
outembimp.em <- scanone(selfC, pheno.col=2, method="imp")

save.image(file="selfC20200708_F.RData")





#!/usr/bin/env Rscript

library(qtl)
setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/selfC")

load("selfC20200708_E.RData")

outepp.em <- scanone(selfC, pheno.col=1)
outemb.em <- scanone(selfC, pheno.col=2)
outeppPermut.em <- scanone(selfC, pheno.col=1, n.perm=1000, verbose=FALSE)
outembPermut.em <- scanone(selfC, pheno.col=2, n.perm=1000, verbose=FALSE)

save.image(file="selfC20200708_G.RData")


#!/usr/bin/env Rscript

library(qtl)
setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/selfC")
selfC <- read.cross("csv","","selfedCforrqtlwithpheno.csv")

selfC <- est.rf(selfC)
selfC <- calc.genoprob(selfC,step=2.5,error.prob=0.01)

save.image(file="selfC20200707.RData")

selfC <- sim.geno(selfC, step=2.5, n.draws=64, error.prob=0.01)
selfC <- sim.geno(selfC, step=2.5, error.prob=0.001, n.draws=256)

selfC_1863 <- sim.geno(subset(selfC, chr=c("Scaffold_1863_HRSCAF_2081")), step=2.5, error.prob=0.001, n.draws=256)


save.image(file="selfC20200707.RData")
