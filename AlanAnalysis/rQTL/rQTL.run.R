### libraries
  library(data.table)
  library(qtl)

### load cross
  AxCF1 <- read.cross("csv","","~/mhp.csv", crosstype="4way", genotypes=NULL)

### reduce to grid
  AxCF1 <- calc.genoprob(AxCF1, stepwidth="fixed")
  AxCF1.sub <- reduce2grid(AxCF1)

  ### marker regression

    mr1 <- scanone(AxCF1, pheno.col=1, method="mr")
    mr2 <- scanone(AxCF1, pheno.col=2, method="mr")
