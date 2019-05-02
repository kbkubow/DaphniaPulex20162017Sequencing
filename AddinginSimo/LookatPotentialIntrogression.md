When Anne Marie constructed a neighbor-joining tree, she found that some DBunk individuals fall out between Dbarb/obstusa and the W and D10 ponds, which was unexpected.
I saw something similar in a PCA. Could this indicate introgression between D. pulex and D. obtusa in DBunk?
We know D. Obtusa is present in DBunk but not really in D8.
First let's record the individuals that fall out in this way in the Tree.

Tree: Spring_2017_DBunk_222, Spring_2017_DBunk_73, December17_2018_D8_1, Spring_2017_DBunk_315, Spring 2017_DBunk_5,
Spring_2017_DBunk_125, Spring_2017_DBunk_21, Spring_2017_DBunk_233, Spring_2017_DBunk_317, Spring_2017_DBunk_242, Spring_2017_DBunk_232,
Spring_2017_DBunk_9, Spring_2017_DBunk_36, Spring_2017_DBunk_3, April_2017_DBunk_185.

Bit disconcerting that the December 2018 individual from D8 is also falling out here. What is going on?

So let's try removing these individuals, then finding fixed differences between D. pulex and D. obtusa.
Then we will go back and see if any of these SNPs are heterozygous in the D. pulex we removed.

Could also try comparing D10/D8 2016/2017 and D. obtusa, and then looking at all individuals in DBunk and seeing if those potentially introgressed individuals look different.
For first pass let's look at just the LD SNPs to make the analysis fo faster. The signal still showed up in the tree based on the LD snps.

Who are the D. obtusa clones? April_2017_Dbarb_11, March20_2018_DBunk_26, March20_2018_DBunk_37, March20_2018_DBunk_42,
March20_2018_DBunk_10, March20_2018_DBunk_18, March20_2018_DBunk_21, March20_2018_DBunk_22, March20_2018_DBunk_23, 
March20_2018_DBunk_38, March20_2018_DBunk_40, March20_2018_DBunk_41, March20_2018_DBunk_43

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
   	genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

#Load individual file
    load("subsetsingleindperSCwcount_20190501.Rdata")
    subsetsingleindperSCwcountids <- subsetsingleindperSCwcount$clone
    temp <- unlist(strsplit(subsetsingleindperSCwcount$clone, split="_"))
		mat <- matrix(temp, ncol=4, byrow=TRUE)
		matdat <- as.data.table(mat)
		subsetsingleindperSCwcount$population <- matdat$V3
		subsetsingleindperSCwcount$year <- matdat$V2
		subsetsingleindperSCwcount$season <- matdat$V1
		subsetsingleindperSCwcount$season2 <- ifelse(tab$season=="April", "Spring", tab$season)
    
#Load LD pruned SNP set
     load("finalsetsnpset01wSimo_20190430.Rdata")

#Pull out D10, D8 and obtusa individuals
     D8D10 <- subsetsingleindperSCwcount[population=="D10" | population=="D8"]
     D8D10no2018 <- D8D10[year!="2018"]
     D8D10no2018$species <- c("pulex")
     Obtusa <- c("April_2017_Dbarb_11", "March20_2018_DBunk_26", "March20_2018_DBunk_37", "March20_2018_DBunk_42",
      "March20_2018_DBunk_10", "March20_2018_DBunk_18", "March20_2018_DBunk_21", "March20_2018_DBunk_22", "March20_2018_DBunk_23", 
      "March20_2018_DBunk_38", "March20_2018_DBunk_40", "March20_2018_DBunk_41", "March20_2018_DBunk_43")
     Obtusadt <- as.data.table(Obtusa)
     colnames(Obtusadt) <- c("clone")
     setkey(Obtusadt, clone)
     setkey(subsetsingleindperSCwcount, clone)
     mObtusa <- merge(Obtusadt, subsetsingleindperSCwcount)
     mObtusa$species <- c("obtusa")
     D8D10no18Obtusa <- rbind(D8D10no2018, mObtusa)
     D8D10no18Obtusaids <- D8D10no18Obtusa$clone
# Set sequence filters
    seqSetFilter(genofile, sample.id=D8D10no18Obtusaids)
    seqSetFilter(genofile, variant.id=finalsetsnpset01)
    
# Call genotypes
    snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		  chr = seqGetData(genofile, "chromosome"),
		  pos = seqGetData(genofile, "position"),
		  dp = seqGetData(genofile, "annotation/info/DP"))

	  het <- t(seqGetData(genofile, "$dosage"))
  	het <- as.data.table(het)
			
	  colnames(het) <- c(seqGetData(genofile, "sample.id"))
	  het$variant.ids <- seqGetData(genofile, "variant.id")
			
	  setkey(het, variant.ids)
	  setkey(snps, variant.ids)
			
	  mhet <- merge(snps, het)

	  sample.ids <- seqGetData(genofile, "sample.id")

#Make long format
	  mhetlong <- melt(mhet, measure.vars=sample.ids, variable.name="clone", value.name="dosage")
	  setkey(mhetlong, clone)
	  setkey(D8D10no18Obtusa, clone)
	  mmhetlong <- merge(D8D10no18Obtusa, mhetlong)
	
#Count instances of variant.ids/SC/dosage	
	dosagecounts <- mmhetlong[, .N, by=list(variant.ids, dosage, species)]

#Remove NAs
	dosagecounts <- dosagecounts[dosage!="NA"]

#Transform to wide format
	doscountwide <- dcast(dosagecounts, variant.ids + species ~ dosage, value.var="N")
	colnames(doscountwide) <- c("variant.ids", "species", "dos0", "dos1", "dos2")
  doscountwide[is.na(dos0),dos0:=0]
	doscountwide[is.na(dos1),dos1:=0]
	doscountwide[is.na(dos2),dos2:=0]
	yetwider <- dcast(doscountwide, variant.ids~species, value.var=c("dos0", "dos1", "dos2"))
  yetwider$obtusa_total <- yetwider$dos0_obtusa + yetwider$dos1_obtusa + yetwider$dos2_obtusa
  yetwider$pulex_total <- yetwider$dos0_pulex + yetwider$dos1_pulex + yetwider$dos2_pulex
  max(yetwider$obtusa_total, na.rm=TRUE) #7
  max(yetwider$pulex_total, na.rm=TRUE) #23
  yetwiderlowmiss <- yetwider[obtusa_total > 5 & pulex_total > 21]
  yetwiderlowmiss$fixed <- ifelse(yetwiderlowmiss$dos0_obtusa/yetwiderlowmiss$obtusa_total==1 &
    yetwiderlowmiss$dos2_pulex/yetwiderlowmiss$pulex_total==1, 1, 0)
  ObtusaPulexfixedids <- yetwiderlowmiss$variant.ids[yetwiderlowmiss$fixed==1]
  
#Now pull out all DBunk and 2018 individuals
  DBunk <- subsetsingleindperSCwcount[population=="DBunk" & year!="2018"]
  all2018 <- subsetsingleindperSCwcount[year=="2018"]
  DBunkand2018 <- rbind(DBunk, all2018)
  DBunkand2018ids <- DBunkand2018$clone

#Set sequence filters
    seqSetFilter(genofile, sample.id=DBunkand2018ids)
    seqSetFilter(genofile, variant.id=ObtusaPulexfixedids)

# Call genotypes
    snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		  chr = seqGetData(genofile, "chromosome"),
		  pos = seqGetData(genofile, "position"),
		  dp = seqGetData(genofile, "annotation/info/DP"))

	  het <- t(seqGetData(genofile, "$dosage"))
  	het <- as.data.table(het)
			
	  colnames(het) <- c(seqGetData(genofile, "sample.id"))
	  het$variant.ids <- seqGetData(genofile, "variant.id")
			
	  setkey(het, variant.ids)
	  setkey(snps, variant.ids)
			
	  mhet <- merge(snps, het)

	  sample.ids <- seqGetData(genofile, "sample.id")

#Make long format
	  mhetlong <- melt(mhet, measure.vars=sample.ids, variable.name="clone", value.name="dosage")
	  setkey(mhetlong, clone)
	  setkey(DBunkand2018, clone)
	  mmhetlong <- merge(DBunkand2018, mhetlong)

#Count instances of clone/dosage	
	  dosagecounts <- mmhetlong[, .N, by=list(clone, dosage)]

#Remove NAs
	  dosagecounts <- dosagecounts[dosage!="NA"]

#Transform to wide format
	  doscountwide <- dcast(dosagecounts, clone ~ dosage, value.var="N")
	  colnames(doscountwide) <- c("clone", "dos0", "dos1", "dos2")
    doscountwide[is.na(dos0),dos0:=0]
	  doscountwide[is.na(dos1),dos1:=0]
	  doscountwide[is.na(dos2),dos2:=0]

# Forgot the obtusa individuals are still in here
    doscountwidenoObtusa <- doscountwide[dos0 < 50000]

```
Not sure this is working quite right. Now let's try just D10 and W individuals versus Obtusa
```

#Pull out D10, W and obtusa individuals
     WD10 <- subsetsingleindperSCwcount[population=="D10" | population=="W1" | population=="W6"]
     WD10$species <- c("pulex")
     Obtusa <- c("April_2017_Dbarb_11", "March20_2018_DBunk_26", "March20_2018_DBunk_37", "March20_2018_DBunk_42",
      "March20_2018_DBunk_10", "March20_2018_DBunk_18", "March20_2018_DBunk_21", "March20_2018_DBunk_22", "March20_2018_DBunk_23", 
      "March20_2018_DBunk_38", "March20_2018_DBunk_40", "March20_2018_DBunk_41", "March20_2018_DBunk_43")
     Obtusadt <- as.data.table(Obtusa)
     colnames(Obtusadt) <- c("clone")
     setkey(Obtusadt, clone)
     setkey(subsetsingleindperSCwcount, clone)
     mObtusa <- merge(Obtusadt, subsetsingleindperSCwcount)
     mObtusa$species <- c("obtusa")
     D10Wno18Obtusa <- rbind(WD10, mObtusa)
     D10Wno18Obtusaids <- D10Wno18Obtusa$clone
# Set sequence filters
    seqSetFilter(genofile, sample.id=D10Wno18Obtusaids)
    seqSetFilter(genofile, variant.id=finalsetsnpset01)
    
# Call genotypes
    snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		  chr = seqGetData(genofile, "chromosome"),
		  pos = seqGetData(genofile, "position"),
		  dp = seqGetData(genofile, "annotation/info/DP"))

	  het <- t(seqGetData(genofile, "$dosage"))
  	het <- as.data.table(het)
			
	  colnames(het) <- c(seqGetData(genofile, "sample.id"))
	  het$variant.ids <- seqGetData(genofile, "variant.id")
			
	  setkey(het, variant.ids)
	  setkey(snps, variant.ids)
			
	  mhet <- merge(snps, het)

	  sample.ids <- seqGetData(genofile, "sample.id")

#Make long format
	  mhetlong <- melt(mhet, measure.vars=sample.ids, variable.name="clone", value.name="dosage")
	  setkey(mhetlong, clone)
	  setkey(D10Wno18Obtusa, clone)
	  mmhetlong <- merge(D10Wno18Obtusa, mhetlong)
	
#Count instances of variant.ids/SC/dosage	
	dosagecounts <- mmhetlong[, .N, by=list(variant.ids, dosage, species)]

#Remove NAs
	dosagecounts <- dosagecounts[dosage!="NA"]

#Transform to wide format
	doscountwide <- dcast(dosagecounts, variant.ids + species ~ dosage, value.var="N")
	colnames(doscountwide) <- c("variant.ids", "species", "dos0", "dos1", "dos2")
  doscountwide[is.na(dos0),dos0:=0]
	doscountwide[is.na(dos1),dos1:=0]
	doscountwide[is.na(dos2),dos2:=0]
	yetwider <- dcast(doscountwide, variant.ids~species, value.var=c("dos0", "dos1", "dos2"))
  yetwider$obtusa_total <- yetwider$dos0_obtusa + yetwider$dos1_obtusa + yetwider$dos2_obtusa
  yetwider$pulex_total <- yetwider$dos0_pulex + yetwider$dos1_pulex + yetwider$dos2_pulex
  max(yetwider$obtusa_total, na.rm=TRUE) #7
  max(yetwider$pulex_total, na.rm=TRUE) #6
  yetwiderlowmiss <- yetwider[obtusa_total > 5 & pulex_total > 4]
  yetwiderlowmiss$fixed <- ifelse(yetwiderlowmiss$dos0_obtusa/yetwiderlowmiss$obtusa_total==1 &
    yetwiderlowmiss$dos2_pulex/yetwiderlowmiss$pulex_total==1, 1, 0)
  ObtusaPulexfixedidsD10W <- yetwiderlowmiss$variant.ids[yetwiderlowmiss$fixed==1]

#Now pull out all DBunk and 2018 individuals
  DBunk <- subsetsingleindperSCwcount[population=="DBunk" & year!="2018"]
  all2018 <- subsetsingleindperSCwcount[year=="2018"]
  DBunkand2018 <- rbind(DBunk, all2018)
  DBunkand2018ids <- DBunkand2018$clone

#Set sequence filters
    seqSetFilter(genofile, sample.id=DBunkand2018ids)
    seqSetFilter(genofile, variant.id=ObtusaPulexfixedidsD10W)

# Call genotypes
    snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		  chr = seqGetData(genofile, "chromosome"),
		  pos = seqGetData(genofile, "position"),
		  dp = seqGetData(genofile, "annotation/info/DP"))

	  het <- t(seqGetData(genofile, "$dosage"))
  	het <- as.data.table(het)
			
	  colnames(het) <- c(seqGetData(genofile, "sample.id"))
	  het$variant.ids <- seqGetData(genofile, "variant.id")
			
	  setkey(het, variant.ids)
	  setkey(snps, variant.ids)
			
	  mhet <- merge(snps, het)

	  sample.ids <- seqGetData(genofile, "sample.id")

#Make long format
	  mhetlong <- melt(mhet, measure.vars=sample.ids, variable.name="clone", value.name="dosage")
	  setkey(mhetlong, clone)
	  setkey(DBunkand2018, clone)
	  mmhetlong <- merge(DBunkand2018, mhetlong)

#Count instances of clone/dosage	
	  dosagecounts <- mmhetlong[, .N, by=list(clone, dosage)]

#Remove NAs
	  dosagecounts <- dosagecounts[dosage!="NA"]

#Transform to wide format
	  doscountwide <- dcast(dosagecounts, clone ~ dosage, value.var="N")
	  colnames(doscountwide) <- c("clone", "dos0", "dos1", "dos2")
    doscountwide[is.na(dos0),dos0:=0]
	  doscountwide[is.na(dos1),dos1:=0]
	  doscountwide[is.na(dos2),dos2:=0]

# Forgot the obtusa individuals are still in here
    doscountwidenoObtusa <- doscountwide[dos0 < 50000]
```
Coming out similar, in that there is one clear outlier, March20_2018_DBunk_5, but others all look pretty similar.
