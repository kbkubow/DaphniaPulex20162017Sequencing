ijob -A berglandlab -p standard -t 05:00:00 -N 1 -n 1 --mem 200000
module load  gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3


library(qtl)
setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/AxCF1/reconstruct")
AxCF1 <- read.cross("csv","","AxCF1genoandphenoreconstruct_sub3.csv", genotypes=NULL, crosstype="4way")
4-way cross

 No. individuals:    23

 No. phenotypes:     4
 Percent phenotyped: 100 100 100 100

 No. chromosomes:    12
     Autosomes:      Scaffold_1863_HRSCAF_2081 Scaffold_1931_HRSCAF_2197
                     Scaffold_2158_HRSCAF_2565 Scaffold_2217_HRSCAF_2652
                     Scaffold_2373_HRSCAF_2879 Scaffold_6786_HRSCAF_7541
                     Scaffold_7757_HRSCAF_8726 Scaffold_9197_HRSCAF_10753
                     Scaffold_9198_HRSCAF_10754 Scaffold_9199_HRSCAF_10755
                     Scaffold_9200_HRSCAF_10757 Scaffold_9201_HRSCAF_10758

 Total markers:      3399
 No. markers:        297 210 264 393 228 345 345 219 351 279 291 177
 Percent genotyped:  100
 Genotypes (%):      AC:24.9  BC:22.9  AD:25.9  BD:26.3  AC/AD:0.0  BC/BD:0.0
                     AC/BC:0.0  AD/BD:0.0  AC/BD:0.0  BC/AD:0.0  not AC:0.0
                     not BC:0.0  not AD:0.0  not BD:0.0

  cg <- comparegeno(AxCF1)
  hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
  rug(cg[lower.tri(cg)])
  wh <- which(cg > 0.9, arr=TRUE)
  wh <- wh[wh[,1] < wh[,2],]
  wh


  #Dropping markers with segregation distortion.
  gt <- geno.table(AxCF1)
  gt[gt$P.value < 0.05/totmar(AxCF1),]
  todrop <- rownames(gt[gt$P.value < 1e-10,])


  AxCF1 <- est.rf(AxCF1)
  plotMap(AxCF1)
  plotRF(AxCF1)

  summaryMap(AxCF1)

  AxCF1 <- calc.genoprob(AxCF1,step=2.5,error.prob=0.01)

  save.image(file="AxCF1_sub3_20200715.RData")

  outemb.em <- scanone(AxCF1, pheno.col=2)
  outemb.hk <- scanone(AxCF1, pheno.col=2, method="hk")
  outembPermut.em <- scanone(AxCF1, pheno.col=2, n.perm=1000, verbose=FALSE)
  outembPermut.hk <- scanone(AxCF1, pheno.col=2, n.perm=100, verbose=FALSE)

  #!/usr/bin/env Rscript

  library(qtl)
  setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/AxCF1/reconstruct")

  load("AxCF1_sub3_20200715.RData")

  outemb.em <- scanone(AxCF1, pheno.col=2)
  outembPermut.em <- scanone(AxCF1, pheno.col=2, n.perm=100, verbose=FALSE)

  save.image(file="AxCF1_sub3_20200715_embscaneone.RData")


  #!/usr/bin/env Rscript

  library(qtl)
  setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/AxCF1/reconstruct")

  load("AxCF1_sub3_20200715.RData")

  outepp.em <- scanone(AxCF1, pheno.col=1)
  outeppPermut.em <- scanone(AxCF1, pheno.col=1, n.perm=100, verbose=FALSE)

  save.image(file="AxCF1_sub3_20200715_eppscaneone.RData")



  #!/usr/bin/env Rscript

  library(qtl)
  setwd("/scratch/kbb7sh/Daphnia/MappingDecember2019/OneLiters/AxCF1/reconstruct")

  load("AxCF1_sub3_20200715.RData")

  AxCF1 <- sim.geno(AxCF1, step=2.5, n.draws=64, error.prob=0.01)


  save.image(file="AxCF1_sub3_20200715_B.RData")
