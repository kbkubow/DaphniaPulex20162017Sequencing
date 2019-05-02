Let's make a neighbor joining tree using LD snps and parsed set of superclones such that each superclone is only represented once regardless of season/year/pond.
```
# Load libraries
    library(vcfR)
    library(ape)
    library(gdsfmt) 
    library(SNPRelate)
    library(data.table)
    library(ggplot2)
    library(foreach)
    library(lattice)
    library(tidyr)
    library(SeqArray)
    library(tidyverse)
    
#Open genotype file
    genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

  #Load parsed superclone file
    load("subsetsingleindperSCwcount_20190501.Rdata")
    subsetsingleindperSCwcountids <- subsetsingleindperSCwcount$clone
   
 #Load SNP file
    load("finalsetsnpset01wSimo_20190430.Rdata")
    
 #Set sequence filters
    seqSetFilter(genofile, sample.id=subsetsingleindperSCwcountids)
    seqSetFilter(genofile, variant.id=finalsetsnpset01)
    
 #Output new VCF
    seqGDS2VCF(genofile, "20161718LDpruned_20190501.vcf")

 # Read VCF into vcfR
    filtvcf <- read.vcfR("20161718LDpruned_20190501.vcf", cols=NULL, checkFile=TRUE, verbose=TRUE)

 #Convert to DNABIN
    LDpruned20161718_20190501DNABIN <- vcfR2DNAbin(filtvcf, 
          extract.indels=TRUE, consensus=TRUE, extract.haps=FALSE, verbose=TRUE)
          
    save(LDpruned20161718_20190501DNABIN, file="LDpruned20161718_20190501DNABIN.Rdata")
    
 # Calculate distance matrix in ape
    LDpruned20161718_20190501DNABINdist <- dist.gene(LDpruned20161718_20190501DNABIN, method="pairwise")

    save(LDpruned20161718_20190501DNABINdist, file="LDpruned20161718_20190501DNABINdist_20190105.Rdata")
    
 # Make neighbor-joining tree

    LDpruned20161718_20190501DNABINdistnjs <- njs(LDpruned20161718_20190501DNABINdist)

    save(LDpruned20161718_20190501DNABINdistnjs, file="LDpruned20161718_20190501DNABINdistnjs_20190501.Rdata")

    plot(root(LDpruned20161718_20190501DNABINdistnjs, 39))
