Figuring out how to calculate LD

```
#!/usr/bin/env Rscript

### libraries
        library(gdsfmt) 
        library(SNPRelate)
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(SeqArray)
        library(tidyverse)

### open geno file
		genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### open subset file
    load("subsetindwithcountsctotal_20190501.Rdata")
    
### Pull out pond/year/etc info and get list of IDs for 2017, DBunk, April clones
    temp <- unlist(strsplit(subsetindwithcountsc$clone, split="_"))
		mat <- matrix(temp, ncol=4, byrow=TRUE)
		matdat <- as.data.table(mat)
		subsetindwithcountsc$population <- matdat$V3
		subsetindwithcountsc$year <- matdat$V2
		subsetindwithcountsc$season <- matdat$V1
    
    DBunk2017April <- subsetindwithcountsc[population=="DBunk" & year=="2017" & season!="May"]
    
    DBunk2017AprilIDs <- DBunk2017April$clone
    
    seqSetFilter(genofile, sample.id=DBunk2017AprilIDs)

    snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
			chr = seqGetData(genofile, "chromosome"),
			pos = seqGetData(genofile, "position"),
			dp = seqGetData(genofile, "annotation/info/DP"))

		
### Now load LD pruned SNPs
    load("finalsetsnpset01dtwSimo_20190430.Rdata")
    
    finalsetsnpset01 <- finalsetsnpset01dt$finalsetsnpset01

    seqSetFilter(genofile, variant.id=finalsetsnpset01)

   snpsLD <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
			chr = seqGetData(genofile, "chromosome"),
			pos = seqGetData(genofile, "position"),
			dp = seqGetData(genofile, "annotation/info/DP"))


### Now initial attempt at calculating LD in windows
   
    load("finalsnpstousewSimoids_20190430.Rdata")

    LDslide20 <- snpgdsLDMat(genofile, sample.id=DBunk2017AprilIDs, snp.id=finalsetsnpset01)
    
    image(t(LDslide20$LD^2), col=terrain.colors(16))

    LDmat <- LDslide20$LD

        ### make the LD matrix long form
                LDlong <- as.data.table(melt(LDmat))


### Ignore below this line.
    LDslide20tmp <- snpgdsLDMat(genofile, sample.id=DBunk2017AprilIDs, snp.id=tmpids, verbose=TRUE)
  LDmattmp <- LDslide20tmp$LD

        ### make the LD matrix long form
                LDlongtmp <- as.data.table(melt(LDmattmp))
 ggplot(LDlongtmp[value>0], aes(Var1,Var2, fill=value)) + geom_raster() + scale_fill_viridis()
