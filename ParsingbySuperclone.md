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
		scs$season <- ifelse(scs$season=="April", "Spring", scs$season)
    
  #Remove two male libraries and SM libraries
    scs <- scs[season!="Lab" & clone!="May_2017_D8_770SM" & clone!="April_2017_D8_515R" & clone!="May_2017_D8_773SM" &
      clone!="Spring_2017_DBunk_116SM" & clone!="Spring_2017_DBunk_347SM" & clone!="May_2017_D8_731SM"]
  
  #Make variable for superclone/year/pond
    scs$supercloneyearpond <- paste(scs$SC, scs$year, scs$population, sep="")
    
  #Based on this new variable, make a list of individuals where a random representative is chosen for each superclone/year/pond
    m_nouniqc <- scs[SC!="OO"]
		m_uniqc <- scs[SC=="OO"]

		m_nouniqc$count <- ifelse(m_nouniqc$year > 0, 1, 0)
		m_uniqc$count <- ifelse(m_uniqc$year > 0, 1, 0)

		m_nouniqc <- m_nouniqc[with(m_nouniqc, order(supercloneyearpond)), ]
		
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
		
