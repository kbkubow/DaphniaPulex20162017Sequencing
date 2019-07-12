Using the parsed set of individuals and LD pruned filtered SNPs, let's try doing some PCA analyses. First let's try with single representatives of each superclone.
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
    
#Load LD pruned SNP set
     load("finalsetsnpset01wSimo_20190430.Rdata")
     

#Set parameters
    ### set some global parameters
		  maf <- 0.01
		  missing.rate <- 0.15
		  threads <- 10

		#### Running pca

		pca <- snpgdsPCA(genofile, snp.id=finalsetsnpset01, sample.id=subsetsingleindperSCwcountids, autosome.only=FALSE, maf=maf, 
									missing.rate=missing.rate, num.thread=threads)
                  #Working space: 132 samples, 134,629 SNVs

		pc.percent <- pca$varprop*100
		head(round(pc.percent, 2))
    #[1] 72.74  5.86  1.74  1.64  0.90  0.75

	
		tab <- data.frame(clone = pca$sample.id,
			EV1 = pca$eigenvect[,1],
			EV2 = pca$eigenvect[,2],
			stringsAsFactors = FALSE)
		head(tab)
		
		temp <- unlist(strsplit(tab$clone, split="_"))
		mat <- matrix(temp, ncol=4, byrow=TRUE)
		matdat <- as.data.table(mat)
		tab$population <- matdat$V3
		tab$year <- matdat$V2
		tab$season <- matdat$V1
		tab$popseason <- paste(tab$population,"_",tab$season, sep="")
		tab$season2 <- ifelse(tab$season=="April", "Spring", tab$season)
		tab$seasonyear <- paste(tab$year,"_",tab$season2, sep="")

		tab <- as.data.table(tab)
		setkey(tab, clone)
		setkey(subsetsingleindperSCwcount, clone)
		mtab <- merge(tab, subsetsingleindperSCwcount)
	
		ggplot(data=mtab, aes(x=EV1, y=EV2, color=popseason)) + geom_point(aes(size=count))
 
 		ggplot(data=mtab[EV1 < 0.2 & EV2 < 0.2], aes(x=EV1, y=EV2, color=popseason)) + geom_point(aes(size=count))
    


    
```

See clear difference between the D. obtusa and D. pulex individuals. Apparently March20_2018_DBunk_10 is also an obtusa, in addition to the superclone in DBarb and DBunk.
Let's try running the PCA again removing all D. obtusa individuals.
```
   #Remove all D. obtusa individuals
    noobstusa <- mtab$clone[mtab$EV1<0.3]
    
   #Rerun PCA
    		pca <- snpgdsPCA(genofile, snp.id=finalsetsnpset01, sample.id=noobstusa, autosome.only=FALSE, maf=maf, 
									missing.rate=missing.rate, num.thread=threads)
        #Working space: 130 samples, 55,803 SNVs

		pc.percent <- pca$varprop*100
		head(round(pc.percent, 2))
	  #[1] 22.24  6.71  6.13  3.35  2.81  2.06

		tab <- data.frame(clone = pca$sample.id,
			EV1 = pca$eigenvect[,1],
			EV2 = pca$eigenvect[,2],
			stringsAsFactors = FALSE)
		head(tab)
		
		temp <- unlist(strsplit(tab$clone, split="_"))
		mat <- matrix(temp, ncol=4, byrow=TRUE)
		matdat <- as.data.table(mat)
		tab$population <- matdat$V3
		tab$year <- matdat$V2
		tab$season <- matdat$V1
		tab$popseason <- paste(tab$population,"_",tab$season, sep="")
		tab$season2 <- ifelse(tab$season=="April", "Spring", tab$season)
		tab$seasonyear <- paste(tab$year,"_",tab$season2, sep="")

		tab <- as.data.table(tab)
		setkey(tab, clone)
		setkey(subsetsingleindperSCwcount, clone)
		mtab <- merge(tab, subsetsingleindperSCwcount)
	
		ggplot(data=mtab, aes(x=EV1, y=EV2, color=popseason)) + geom_point(aes(size=count))
```
![PCAnoObtusa](PCAnoObtusa.tiff)

See clear separation of W1, W6, D10, and D8 and surrounding ponds. Within D10, Fall_2016_D10_41 and the superclone containing Fall_2016_D10_45 and Fall_2016_D10_84 fall out distinct from the superclone containing the remainder of D10 individuals.
Let's try running the PCA on just the D8 and surround ponds (drop W1, W6, and D10).
```
   #Remove W1, W6, and D10 individuals
    noobtusaW1W6D10 <- mtab$clone[mtab$population!="W1" & mtab$population!="W6" & mtab$population!="D10"]
    
       #Rerun PCA
    		pca <- snpgdsPCA(genofile, snp.id=finalsetsnpset01, sample.id=noobtusaW1W6D10, autosome.only=FALSE, maf=maf, 
									missing.rate=missing.rate, num.thread=threads)
        #Working space: 124 samples, 39,475 SNVs

		pc.percent <- pca$varprop*100
		head(round(pc.percent, 2))
	  #[1] 9.99 4.93 3.48 3.36 2.94 2.54

		tab <- data.frame(clone = pca$sample.id,
			EV1 = pca$eigenvect[,1],
			EV2 = pca$eigenvect[,2],
			stringsAsFactors = FALSE)
		head(tab)
		
		temp <- unlist(strsplit(tab$clone, split="_"))
		mat <- matrix(temp, ncol=4, byrow=TRUE)
		matdat <- as.data.table(mat)
		tab$population <- matdat$V3
		tab$year <- matdat$V2
		tab$season <- matdat$V1
		tab$season2 <- ifelse(tab$season=="Spring" & tab$year=="2017", "April", tab$season)
		tab$popseason <- paste(tab$population,"_",tab$season2, sep="")
		tab$seasonyear <- paste(tab$year,"_",tab$season2, sep="")

		tab <- as.data.table(tab)
		setkey(tab, clone)
		setkey(subsetsingleindperSCwcount, clone)
		mtab <- merge(tab, subsetsingleindperSCwcount)
	
		ggplot(data=mtab, aes(x=EV1, y=EV2, color=popseason)) + geom_point(aes(size=count))
		ggplot(data=mtab[population=="D8"], aes(x=EV1, y=EV2, color=popseason)) + geom_point(aes(size=count))
		ggplot(data=mtab[population=="D8" | population=="DBunk"], aes(x=EV1, y=EV2, color=popseason)) + geom_point(aes(size=count))
		ggplot(data=mtab[year!="2012" & year!="2016" &population=="D8"], aes(x=EV1, y=EV2, color=year)) + geom_point(aes(size=count)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
		
```
![PCAnoObtusaW1W6D10](PCAnoObtusaW1W6D10.tiff)
![PCAnoObtusaW1W6D10graphonlyD8](PCAnoObtusaW1W6D10graphonlyD8.tiff)

At least from this PCA, it looks like the March 2018 individuals are to a large extent F1 hybrids of A and B superclones. Get somewhat different patterns if I graph D8 an surrounding individuals but include Obtusa, D10, or W ponds.
Might be worth trying to construct a tree on the same set of individuals with the same LD pruned SNPs.



    
       
