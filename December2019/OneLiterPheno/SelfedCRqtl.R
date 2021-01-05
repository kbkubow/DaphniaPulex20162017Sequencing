ijob -A berglandlab -p standard -t 05:00:00 -N 1 -n 1 --mem 200000
module load  gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3

#!/usr/bin/env Rscript

library(qtl)
setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/selfC")
selfC <- read.cross("csv","","selfedCforrqtlwithpheno_sub.csv")
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

selfC <- switch.order(selfC, chr="Scaffold_1931_HRSCAF_2197", c(1:70, 115:120, 71:114, 121:160), error.prob=0.005)
selfC <- switch.order(selfC, chr="Scaffold_1931_HRSCAF_2197", c(1:50, 53:68, 71:76, 51:52, 69:70, 77:160), error.prob=0.005)
selfC <- switch.order(selfC, chr="Scaffold_1931_HRSCAF_2197", c(1:76, 84, 77:83, 85:160), error.prob=0.005)

rip1 <- ripple(selfC, chr="Scaffold_1931_HRSCAF_2197", window=7)

selfC <- switch.order(selfC, chr="Scaffold_1931_HRSCAF_2197", c(1:65, 67:72, 66, 73:160), error.prob=0.005)

rip1 <- ripple(selfC, chr="Scaffold_1931_HRSCAF_2197", window=7)

# Let's just leave this for now.

plotRF(selfC, chr="Scaffold_1863_HRSCAF_2081")

# Looks fine

plotRF(selfC, chr="Scaffold_2158_HRSCAF_2565")

selfC <- switch.order(selfC, chr="Scaffold_2158_HRSCAF_2565", c(1:90, 102, 91:101, 103:120), error.prob=0.005)

# OK

plotRF(selfC, chr="Scaffold_2217_HRSCAF_2652")

#Leaving for now

plotRF(selfC, chr="Scaffold_2373_HRSCAF_2879")

selfC <- drop.markers (selfC, find.marker(selfC, "Scaffold_2373_HRSCAF_2879", index=1))

#Ok good

plotRF(selfC, chr="Scaffold_6786_HRSCAF_7541")

#Fine

plotRF(selfC, chr="Scaffold_7757_HRSCAF_8726")

selfC <- switch.order(selfC, chr="Scaffold_7757_HRSCAF_8726", c(1:60, 64:93, 101:120, 61:63, 94:100), error.prob=0.005)
selfC <- switch.order(selfC, chr="Scaffold_7757_HRSCAF_8726", c(1:67, 78:87, 68:77, 88:120), error.prob=0.005)






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


outemb.em <- scanone(selfC, pheno.col=2)
outembPermut.em <- scanone(selfC, pheno.col=2, n.perm=1000, verbose=FALSE)

save.image("selfC20200708_I.RData")

summary(outembPermut.em, alpha=c(0.20, 0.05))
LOD thresholds (100 permutations)
     lod
20% 3.65
5%  4.06

summary(outemb.em, perms=outembPermut.em, alpha=0.05)
# None

summary(outemb.em, perms=outembPermut.em, alpha=0.10)
#None

summary(outemb.em, perms=outembPermut.em, alpha=0.15)
#None

summary(outemb.em, perms=outembPermut.em, alpha=0.30)
#chr     pos  lod
#Scaffold_9197_HRSCAF_10753_1097810 Scaffold_9197_HRSCAF_10753 1097810 3.54


outepp.em <- scanone(selfC, pheno.col=1)
outeppPermut.em <- scanone(selfC, pheno.col=1, n.perm=100, verbose=FALSE)

save("selfC20200708_J.RData")

summary(outeppPermut.em, alpha=c(0.20, 0.05))
lod
20% 3.31
5%  3.79

summary(outepp.em, perms=outeppPermut.em, alpha=0.25)
chr     pos  lod
Scaffold_2373_HRSCAF_2879_5741018 Scaffold_2373_HRSCAF_2879 5741018 3.21


#!/usr/bin/env Rscript

library(qtl)
setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/selfC")

load("selfC20200708_I.RData")
outembtwo <- scantwo(selfC, pheno.col=2, verbose=FALSE)

save.image("selfC20200708_I_B.RData")

#!/usr/bin/env Rscript

library(qtl)
setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/selfC")

load("selfC20200708_J.RData")
outepptwo <- scantwo(selfC, pheno.col=1, verbose=FALSE)

save.image("selfC20200708_J_B.RData")
