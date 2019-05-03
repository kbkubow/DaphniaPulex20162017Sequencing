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

        load("finalsnpstousewSimoids_20190430.Rdata")
        
        load("secondsampstokeepwSimo_20190430.Rdata")
        
 # Load superclone file
 
        scs <- fread("Superclones20161718updated_20190501")

# Set sequence filters

        seqSetFilter(genofile, variant.id=finalsnpstousewSimoids)
        seqSetFilter(genofile, sample.id=secondsampstokeep)
        
# Get median read depths

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
           
          write.table(readdepthrow.ag, file="medianreaddepthvcfD_20190503", row.names=FALSE, quote=FALSE, sep="\t")

# Add in superclone info
          setkey(readdepthrow.ag, clone)
          setkey(scs, clone)
          scswrd <- merge(scs, readdepthrow.ag)
          
# Add in year/pop/season info
        temp <- unlist(strsplit(scswrd$clone, split="_"))
        mat <- matrix(temp, ncol=4, byrow=TRUE)
        matdat <- as.data.table(mat)
        scswrd$pop <- matdat$V3
        scswrd$year <- matdat$V2
        scswrd$season <- matdat$V1
        scswrd$seasonB <- ifelse(scswrd$season=="Spring" & scswrd$year=="2017", "April", scswrd$season)
        
# Pull out SCAs from 2017April, 2017May, and 2017March
        SCAD8April2017medrd10 <- scswrd[SC=="A" & season=="April" & year=="2017" & pop=="D8" & medrd > 9]
        SCAD8May2017medrd10 <- scswrd[SC=="A" & season=="May" & year=="2017" & pop=="D8" & medrd > 9]
        SCAD8March2018medrd10 <- scswrd[SC=="A" & season=="March20" & year=="2018" & pop=="D8" & medrd > 9]
        SCAD8March2018all <- scswrd[SC=="A" & season=="March20" & year=="2018" & pop=="D8"]
        totalAsforgeneconversion <- rbind(SCAD8April2017medrd10, SCAD8May2017medrd10, SCAD8March2018all)
        write.table(totalAsforgeneconversion, file="totalAsforgeneconversion_20190503", sep="\t", row.names=FALSE, quote=FALSE)
        
        
