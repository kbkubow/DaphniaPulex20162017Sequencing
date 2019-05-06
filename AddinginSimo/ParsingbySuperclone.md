Let's make a parsed set of individuals such that each superclone is only represented once per year/pond. Also make one where every superclone is only represented once. Let's start with superclone per year/pond.
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

  #Open superclone file of reasonable coverage individuals
    scs <- fread("Superclones20161718updated_20190501")

 #Fix one DCat clone
    scs$clone <- str_replace(scs$clone, "cat", "Cat")
 
  #Add columns to superclone file for pond, year, season.
    temp <- unlist(strsplit(scs$clone, split="_"))
		mat <- matrix(temp, ncol=4, byrow=TRUE)
		matdat <- as.data.table(mat)
		scs$population <- matdat$V3
		scs$year <- matdat$V2
		scs$season <- matdat$V1
		scs$season <- ifelse(scs$season=="Spring", "April", scs$season)
    
  #Remove two male libraries and SM libraries
    scs <- scs[season!="Lab" & clone!="May_2017_D8_770SM" & clone!="April_2017_D8_515R" & clone!="May_2017_D8_773SM" &
      clone!="Spring_2017_DBunk_116SM" & clone!="Spring_2017_DBunk_347SM" & clone!="May_2017_D8_731SM"]
  
  #Make variable for superclone/year/pond/season
    scs$supercloneyearpond <- paste(scs$SC, scs$year, scs$population, scs$season, sep="_")
  
  #Add in median read depths
  	
	load("finalsnpstousewSimoids_20190430.Rdata")
        
        load("secondsampstokeepwSimo_20190430.Rdata")
        
     # Set sequence filters

        seqSetFilter(genofile, sample.id=secondsampstokeep)
        seqSetFilter(genofile, variant.id=finalsnpstousewSimoids)

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
	  
	  readdepthrow.ag <- readdepthrow[,list(medrd=as.double(median(dp1, na.rm=TRUE))), list(clone)]

	setkey(readdepthrow.ag, clone)
	setkey(scs, clone)
	m2 <- merge(readdepthrow.ag, scs)
  
  #Based on this new variable, make a list of individuals where a random representative is chosen for each superclone/year/pond
    		m_nouniqc <- m2[SC!="OO"]
		m_nouniqc <- m_nouniqc[medrd>9]

		m_uniqc <- m2[SC=="OO"]

		m_nouniqc$count <- ifelse(m_nouniqc$year > 0, 1, 0)
		m_uniqc$count <- ifelse(m_uniqc$year > 0, 1, 0)

		m_nouniqc <- m_nouniqc[with(m_nouniqc, order(supercloneyearpond, -medrd)), ]
		
		m_nouniqcountsc <- as.data.table(aggregate(count~supercloneyearpond, m_nouniqc, FUN=sum))

		t.firstc <- m_nouniqc[match(unique(m_nouniqc$supercloneyearpond), m_nouniqc$supercloneyearpond),]
		
		t.firstsubc <- data.table(clone=t.firstc$clone, SC=t.firstc$SC, supercloneyearpond=t.firstc$supercloneyearpond)
		
		setkey(t.firstsubc, supercloneyearpond)
		setkey(m_nouniqcountsc, supercloneyearpond)
		
		mnouniqc <- merge(t.firstsubc, m_nouniqcountsc)
		
		mnouniqBc <- data.table(clone=mnouniqc$clone, SC=mnouniqc$SC, supercloneyearpond=mnouniqc$supercloneyearpond, count=mnouniqc$count)
		
		m_uniqsubc <- data.table(clone=m_uniqc$clone, SC=m_uniqc$SC, supercloneyearpond=m_uniqc$supercloneyearpond, count=m_uniqc$count)
		
		subsetindwithcountsc <- rbind(mnouniqBc, m_uniqsubc)
		
		save(subsetindwithcountsc, file="subsetindwithcountsctotal_20190501.Rdata")
		
    subsetindtouse <- subsetindwithcountsc$clone
    
   
```
This results in a parsed set of 141 individuals. Now let's try choosing just one representative individual per superclone regardless of year, pond, etc.
```
		m_nouniqc <- m_nouniqc[with(m_nouniqc, order(season, -medrd)), ]

		m_nouniqcountsc <- as.data.table(aggregate(count~SC, m_nouniqc, FUN=sum))

		t.firstc <- m_nouniqc[match(unique(m_nouniqc$SC), m_nouniqc$SC),]
		
		t.firstsubc <- data.table(clone=t.firstc$clone, SC=t.firstc$SC, supercloneyearpond=t.firstc$supercloneyearpond)
		
		setkey(t.firstsubc, SC)
		setkey(m_nouniqcountsc, SC)
		
		mnouniqc <- merge(t.firstsubc, m_nouniqcountsc)
		
		mnouniqBc <- data.table(clone=mnouniqc$clone, SC=mnouniqc$SC, supercloneyearpond=mnouniqc$supercloneyearpond, count=mnouniqc$count)
		
		m_uniqsubc <- data.table(clone=m_uniqc$clone, SC=m_uniqc$SC, supercloneyearpond=m_uniqc$supercloneyearpond, count=m_uniqc$count)
		
		subsetindwithcountsc <- rbind(mnouniqBc, m_uniqsubc)
		subsetsingleindperSCwcount <- subsetindwithcountsc

		save(subsetsingleindperSCwcount, file="subsetsingleindperSCwcount_20190501.Rdata")
		
