  Let's make a filtered VCF that contains LD pruned, filtered SNPs for individuals representing each SC only once.
  ```
  
  #Load libraries
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
   	genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

#Load individual file
    load("subsetsingleindperSCwcount_20190425.Rdata")
    subsetsingleindperSCwcountids <- subsetsingleindperSCwcount$clone
    
#Load LD pruned SNP set
     load("finalsetsnpset01_20190423.Rdata")
     
#Set filters
      seqSetFilter(genofile, variant.id=finalsetsnpset01)
      seqSetFilter(genofile, sample.id=subsetsingleindperSCwcountids)
      
#Output vcf with LD pruned, filtered SNPs and for individuals that represent each SC once
      seqGDS2VCF(genofile, "LDprunedfilteredSNPs_eachSConlyone_20190425.vcf")
     
```
Ok, now let's write an Rscript to make a NJS tree.
```
 
#!/usr/bin/env Rscript

### Load libraries

library(vcfR)
library(ape)

### Run program
# Read in VCF

LDprunedSConce <- read.vcfR("LDprunedfilteredSNPs_eachSConlyone_20190425.vcf", cols=NULL, checkFile=TRUE, verbose=TRUE)

# Convert to DNABIN

LDprunedSConceDNABIN <- vcfR2DNAbin(LDprunedSConce, extract.indels=TRUE, consensus=TRUE, extract.haps=FALSE, verbose=TRUE)

save(LDprunedSConceDNABIN, file="LDprunedSConceDNABIN_20190425.Rdata")

# Calculate distance matrix in ape

LDprunedSConceDNABINdnadist <- dist.gene(LDprunedSConceDNABIN, method="pairwise")

save(LDprunedSConceDNABINdnadist, file="LDprunedSConceDNABINdnadist_20190425.Rdata")

# Make neighbor-joining tree

LDprunedSConceDNABINdnadistnjs <- njs(LDprunedSConceDNABINdnadist)

save(LDprunedSConceDNABINdnadistnjs, file="LDprunedSConceDNABINdnadistnjs_20190425.Rdata")

plot(root(LDprunedSConceDNABINdnadistnjs, 38), type="fan")
```


