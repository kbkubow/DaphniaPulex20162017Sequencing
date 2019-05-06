```
#Load libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(foreach)
        library(SeqArray)
        library(doMC)
        registerDoMC(20)
        library(ggplot2)
        library(viridis)
        library(tidyverse)
        
# Load genotype file and snp and sample filter files

        genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
        
        load("finalsetsnpset01wSimo_20190430.Rdata")
        
        load("secondsampstokeepwSimo_20190430.Rdata")

# Run IBS

        ibs <- snpgdsIBS(genofile, snp.id=finalsetsnpset01, sample.id=secondsampstokeep, num.thread=20, maf=0.01,
          missing.rate=0.15, autosome.only = FALSE)
        
        ### a bit of re-formating of the ibs matrix
                ibs.mat <- ibs$ibs
                rownames(ibs.mat) <- ibs$sample.id
                colnames(ibs.mat) <- ibs$sample.id

                ibs.mat[upper.tri(ibs.mat)] <- NA
                
        ### make the IBs matrix long form
                ibs.long <- as.data.table(melt(ibs.mat))
                setnames(ibs.long, names(ibs.long), c("cloneA", "cloneB", "distance"))
                ibs.long <- na.omit(ibs.long)
                
  ### Identify superclones
                ### give temporary labels to super clones based on identity of clone B
                        superclone <- ibs.long[distance>.9865]
                        superclone[,SC.sub := as.numeric(as.factor(cloneB))]

                ### collapse nested superclones
                superclone.o <- foreach(sc.i=unique(superclone$SC.sub), .combine="c")%do%{
                        paste(sort(unique(superclone[cloneB%in%superclone[SC.sub==sc.i]$cloneA]$cloneA)), collapse=";")
                }

                superclone.o.uniq <- unique(superclone.o)

                sc.dt <- foreach(i=superclone.o.uniq, .combine="rbind")%do%{
                        data.table(clone=strsplit(i, ";")[[1]],
                                                superClone.size=length(strsplit(i, ";")[[1]]),
                                                superClone.index=which(i==superclone.o.uniq))
                }
                sc.dt[,superClone.sizeRank := as.numeric(as.factor(rank(-superClone.size, ties="average")))]

                sc.dt <- as.data.table(sc.dt %>% 
                        mutate(SCnum = group_indices_(sc.dt, .dots=c("superClone.sizeRank", "superClone.index"))))

    ### do plot to make sure our head is screwed on correctly
      plot(superClone.sizeRank ~ superClone.size, sc.dt)

    ### label superclones with letters. What do you do when you have more than 26 superclones?
        sc.dt[,SC:=LETTERS[SCnum]]
        sc.dt$SC <- ifelse(sc.dt$SCnum==27, "AA", ifelse(sc.dt$SCnum==28, "AB", ifelse(sc.dt$SCnum==29, "AC",
            ifelse(sc.dt$SCnum==30, "AD", ifelse(sc.dt$SCnum==31, "AE", ifelse(sc.dt$SCnum=="32", "AF", sc.dt$SC))))))

 
    ### rename singletone individuals to "OO" to follow Karen's convention
        sc.dt[superClone.size==1, SC:="OO"]

    ### Adding superclone ids to distance file
        sc.dtsub <- data.table(clone=sc.dt$clone, SC=sc.dt$SC)
        sc.dtsubB <- sc.dtsub
        colnames(sc.dtsub) <- c("cloneA", "SCA")
        setkey(sc.dtsub, cloneA)
        setkey(superclone, cloneA)
        msc <- merge(sc.dtsub, superclone)
        colnames(sc.dtsubB) <- c("cloneB", "SCB")
        setkey(sc.dtsubB, cloneB)
        setkey(msc, cloneB)
        mmsc <- merge(sc.dtsubB, msc)
        
 #Remove entires of clones to themselves
        mscnoident <- mmsc[cloneB!=cloneA]
        
 #Pull out just superclone A individuals
        mscnoidentAs <- mscnoident[SCB=="A" & SCA=="A"]
        
 #Adding pond, year, and season designations for clones A and B
        mscnoidentAs$cloneB <- as.character(mscnoidentAs$cloneB)
        temp <- unlist(strsplit(mscnoidentAs$cloneB, split="_"))
        mat <- matrix(temp, ncol=4, byrow=TRUE)
        matdat <- as.data.table(mat)
        mscnoidentAs$popB <- matdat$V3
        mscnoidentAs$yearB <- matdat$V2
        mscnoidentAs$seasonB <- matdat$V1
        mscnoidentAs$cloneA <- as.character(mscnoidentAs$cloneA)
	temp <- unlist(strsplit(mscnoidentAs$cloneA, split="_"))
        mat <- matrix(temp, ncol=4, byrow=TRUE)
        matdat <- as.data.table(mat)
        mscnoidentAs$popA <- matdat$V3
        mscnoidentAs$yearA <- matdat$V2
        mscnoidentAs$seasonA <- matdat$V1

#Add a variable to pull out year combinations
        mscnoidentAs$yearA_yearB <- paste(mscnoidentAs$yearA,mscnoidentAs$yearB, sep="_")

#Remove single male libraries
        mscnoidentAsnomale <- mscnoidentAs[yearA!="2019" & yearB!="2019"]
        
#Consolidate year designations
        mscnoidentAsnomale$yearA_yearB_B <- ifelse(mscnoidentAsnomale$yearA_yearB=="2017_2016", 
	"2016_2017", ifelse(mscnoidentAsnomale$yearA_yearB=="2018_2017", "2017_2018", 
	mscnoidentAsnomale$yearA_yearB))
		        
#Let's add in median read depths for each sample to see if that is causing some of the outliers
#Set sequence filters
        seqSetFilter(genofile, sample.id=secondsampstokeep)
        seqSetFilter(genofile, variant.id=finalsetsnpset01)

#Pull out median read depths
        snpsG <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
					chr = seqGetData(genofile, "chromosome"),
					pos = seqGetData(genofile, "position"),
					dp = seqGetData(genofile, "annotation/info/DP"))

          dp <- t(seqGetData(genofile, "annotation/format/DP")$data)
          dp <- as.data.table(dp)
          dp[,variant.ids:=snpsG$variant.ids]
	
		      sampleids <- as.data.table(seqGetData(genofile, "sample.id"))
		      variantidhead <- as.data.table("variant.ids")
		      totalids <- rbind(sampleids, variantidhead)
		      colnames(dp) <- c(totalids$V1)
	
	        setkey(dp, variant.ids)
	        setkey(snpsG, variant.ids)
	    
	        m <- merge(snpsG, dp)

		      readdepthrow <- melt(m, measure.vars=sampleids, variable.name="clone", value.name="dp")

		      readdepthrow.ag <- readdepthrow[,list(medrd=as.double(median(dp1, na.rm=TRUE))),
					 list(clone)]

# Add in read depth information to superclone A IBS data table
          colnames(readdepthrow.ag) <- c("cloneA", "medrdA")
          setkey(readdepthrow.ag, cloneA)
          setkey(mscnoidentAsnomale, cloneA)
          m <- merge(readdepthrow.ag, mscnoidentAsnomale)
          colnames(readdepthrow.ag) <- c("cloneB", "medrdB")
          setkey(readdepthrow.ag, cloneB)
          setkey(m, cloneB)
          m2 <- merge(readdepthrow.ag, m)
          m2$avgrd <- (m2$medrdB+m2$medrdA)/2

#Look at the relationship between average read depth and distance
          ggplot(data=m2, aes(x=avgrd, y=distance)) + geom_point()

#Remove 2017_May_D8_770SM which seems to be an outlier
	 m2no770SM <- m2[cloneA!="2017_May_D8_770SM" & cloneB!="2017_May_D8_770SM"]

#There is a relationship. It seems like low read depth individuals are causing some of the IBS outliers.
          ggplot(data=m2, aes(x=yearA_yearB_B, y=distance)) + geom_violin() + ylim(0.985,1)
          ggplot(data=m2[avgrd>15], aes(x=yearA_yearB_B, y=distance)) + geom_violin()
          
#But perhaps individual median read depth is better to control for than average read depth
          ggplot(data=m2[medrdA>8 & medrdB>8], aes(x=yearA_yearB_B, y=distance)) + geom_violin() + ylim(0.985,1)
          ggplot(data=m2[medrdA>9 & medrdB>9], aes(x=yearA_yearB_B, y=distance)) + geom_violin() + ylim(0.985,1)
          ggplot(data=m2[medrdA>14 & medrdB>14], aes(x=yearA_yearB_B, y=distance)) + geom_violin() + ylim(0.985,1)

	m2$seasonA_B <- ifelse(m2$seasonA=="Spring" & m2$yearA=="2017", "April", m2$seasonA)
	m2$seasonB_B <- ifelse(m2$seasonB=="Spring" & m2$yearB=="2017", "April", m2$seasonB)
	m2$seasonA_yearA <- paste(m2$seasonA_B, m2$yearA, sep="_")
	m2$seasonB_yearB <- paste(m2$seasonB_B, m2$yearB, sep="_")
	m2$seasonyearcomp <- paste(m2$seasonA_yearA, m2$seasonB_yearB, sep="_")

#Consolidate seasonyearcomp
	m2$seasonyearcompB <- ifelse(m2$seasonyearcomp=="April_2017_Spring_2016", "Spring_2016_April_2017",
		ifelse(m2$seasonyearcomp=="March20_2018_April_2017", "April_2017_March20_2018",
		ifelse(m2$seasonyearcomp=="May_2017_April_2017", "April_2017_May_2017",
		m2$seasonyearcomp)))
	m2$seasonyearcompB <- factor(m2$seasonyearcompB, levels = m2$seasonyearcompB[order(m2$yearA, m2$yearB, m2$seasonA_B)])
       
        ggplot(data=m2, aes(x=seasonyearcompB, y=distance)) + geom_violin() + ylim(0.985,1)
        ggplot(data=m2[medrdA>9 & medrdB>9], aes(x=seasonyearcompB, y=distance)) + geom_violin() + ylim(0.985,1)
        ggplot(data=m2[medrdA>14 & medrdB>14], aes(x=seasonyearcompB, y=distance)) + geom_violin() + ylim(0.985,1)
 		
	
	
###What happens if we look at superclone B as comparison? Don't have 2018, so not as useful. But let's just take a look.

#Pull out just superclone A individuals
        mscnoidentBs <- mscnoident[SCB=="B" & SCA=="B"]
        
 #Adding pond, year, and season designations for clones A and B
        mscnoidentBs$cloneB <- as.character(mscnoidentBs$cloneB)
        temp <- unlist(strsplit(mscnoidentBs$cloneB, split="_"))
        mat <- matrix(temp, ncol=4, byrow=TRUE)
        matdat <- as.data.table(mat)
        mscnoidentBs$popB <- matdat$V3
        mscnoidentBs$yearB <- matdat$V2
        mscnoidentBs$seasonB <- matdat$V1
        mscnoidentBs$cloneA <- as.character(mscnoidentBs$cloneA)
	temp <- unlist(strsplit(mscnoidentBs$cloneA, split="_"))
        mat <- matrix(temp, ncol=4, byrow=TRUE)
        matdat <- as.data.table(mat)
        mscnoidentBs$popA <- matdat$V3
        mscnoidentBs$yearA <- matdat$V2
        mscnoidentBs$seasonA <- matdat$V1

#Add a variable to pull out year combinations
        mscnoidentBs$yearA_yearB <- paste(mscnoidentBs$yearA,mscnoidentBs$yearB, sep="_")

#Remove single male libraries
        mscnoidentBsnomales <- mscnoidentBs[yearA!="2019" & yearB!="2019"]
        
#Consolidate year designations
        mscnoidentBsnomales$yearA_yearB_B <- ifelse(mscnoidentBsnomales$yearA_yearB=="2017_2016", 
		"2016_2017", mscnoidentBsnomales$yearA_yearB)
		        
#Let's add in median read depths for each sample to see if that is causing some of the outliers

# Add in read depth information to superclone A IBS data table
          colnames(readdepthrow.ag) <- c("cloneA", "medrdA")
          setkey(readdepthrow.ag, cloneA)
          setkey(mscnoidentBsnomales, cloneA)
          mB <- merge(readdepthrow.ag, mscnoidentBsnomales)
          colnames(readdepthrow.ag) <- c("cloneB", "medrdB")
          setkey(readdepthrow.ag, cloneB)
          setkey(mB, cloneB)
          m2B <- merge(readdepthrow.ag, mB)
          m2B$avgrd <- (m2B$medrdB+m2B$medrdA)/2

#Look at the relationship between average read depth and distance
          ggplot(data=m2B, aes(x=avgrd, y=distance)) + geom_point()

#Remove 2017_May_D8_773SM which seems to be an outlier
	 m2Bno773SM <- m2B[cloneA!="May_2017_D8_773SM" & cloneB!="May_2017_D8_773SM"]

#There is a relationship. It seems like low read depth individuals are causing some of the IBS outliers.
          ggplot(data=m2B, aes(x=yearA_yearB_B, y=distance)) + geom_violin() + ylim(0.985,1)
          ggplot(data=m2B[avgrd>15], aes(x=yearA_yearB_B, y=distance)) + geom_violin()
          
#But perhaps individual median read depth is better to control for than average read depth
          ggplot(data=m2B[medrdA>8 & medrdB>8], aes(x=yearA_yearB_B, y=distance)) + geom_violin() + ylim(0.985,1)
          ggplot(data=m2B[medrdA>9 & medrdB>9], aes(x=yearA_yearB_B, y=distance)) + geom_violin() + ylim(0.985,1)
          ggplot(data=m2B[medrdA>14 & medrdB>14], aes(x=yearA_yearB_B, y=distance)) + geom_violin() + ylim(0.985,1)

```
