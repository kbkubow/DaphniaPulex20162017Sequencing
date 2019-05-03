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
        
#Remove multiple comparisons among years
        mscnoidentAsnomalesing <- mscnoidentAsnomale[yearA_yearB!="2017_2016" & yearA_yearB!="2018_2016" & yearA_yearB!="2018_2017"]
       
#Remove D8_770SM library, which seems to be an outlier
        mscnoidentAsnomalesingno770SM <0 mscnoidentAsnomalesing[cloneA!="May_2017_D8_770SM" & cloneB!="May_2017_D8_770SM"]
        
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
          setkey(mscnoidentAsnomalesingno770SM, cloneA)
          m <- merge(readdepthrow.ag, mscnoidentAsnomalesingno770SM)
          colnames(readdepthrow.ag) <- c("cloneB", "medrdB")
          setkey(readdepthrow.ag, cloneB)
          setkey(m, cloneB)
          m2 <- merge(readdepthrow.ag, m)
          m2$avgrd <- (m2$medrdB+m2$medrdA)/2
 
#Look at the relationship between average read depth and distance
          ggplot(data=m2, aes(x=avgrd, y=distance)) + geom_point()
          
#There is a relationship. It seems like low read depth individuals are causing some of the IBS outliers.
          ggplot(data=m2, aes(x=yearA_yearB, y=distance)) + geom_violin()
          ggplot(data=m2[avgrd>15], aes(x=yearA_yearB, y=distance)) + geom_violin()
          
#But perhaps individual median read depth is better to control for than average read depth
          ggplot(data=m2[medrdA>9 & medrdB>9], aes(x=yearA_yearB, y=distance)) + geom_violin()
          ggplot(data=m2[medrdA>14 & medrdB>14], aes(x=yearA_yearB, y=distance)) + geom_violin()
        
