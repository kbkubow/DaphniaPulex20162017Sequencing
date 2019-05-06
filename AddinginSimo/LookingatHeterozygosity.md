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

        load("finalsnpstousewSimoids_20190430.Rdata")
        
        load("secondsampstokeepwSimo_20190430.Rdata")
        
# Set sequence filters

        seqSetFilter(genofile, sample.id=secondsampstokeepwSimo)
        seqSetFilter(genofile, variant.id=finalsnpstousewSimoids)
        
# Load superclone file
      
        scs <- fread("Superclones20161718updated_20190501")
        
# Pull out all A clones

        scsA <- scs[SC=="A"]

#Remove males and SM libraries

        scsAnomaleSM <- scsA[clone!="Lab_2019_D8_349Male" & clone!="May_2017_D8_770SM"]
        scsAnomaleSMids <- scsAnomaleSM$clone
        
# Reset sequence filter

        seqSetFilter(genofile, sample.id=scsAnomaleSMids)
        
# Pull out genotypes

        snpsG <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
		      chr = seqGetData(genofile, "chromosome"),
		      pos = seqGetData(genofile, "position"),
		      dp = seqGetData(genofile, "annotation/info/DP"))

	      het <- t(seqGetData(genofile, "$dosage"))
	      het <- as.data.table(het)
			
	      colnames(het) <- c(seqGetData(genofile, "sample.id"))
	      het$variant.ids <- seqGetData(genofile, "variant.id")
			
	      setkey(het, variant.ids)
	      setkey(snpsG, variant.ids)
			
	      mhet <- merge(snpsG, het)

	      mhetlong <- melt(mhet, measure.vars=scsAnomaleSMids, variable.name="clone", value.name="dosage")

#Count instances of variant.ids/clone/dosage	
	      genoclone <- mhetlong[, .N, by=list(clone, dosage)]

#Remove NAs
	genoclone <- genoclone[dosage!="NA"]

#Transform to wide format
	genoclonewide <- dcast(genoclone, clone ~ dosage, value.var="N")
	colnames(genoclonewide) <- c("clone", "dos0", "dos1", "dos2")
	genoclonewide$total <- genoclonewide$dos1 + genoclonewide$dos2 + genoclonewide$dos0
	genoclonewide$prophet <- genoclonewide$dos1/genoclonewide$total
	genoclonewide$hetoverhomalt <- genoclonewide$dos1/(genoclonewide$dos0+genoclonewide$dos1)

#Add in median read depth

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
	  setkey(genoclonewide, clone)
	  m <- merge(readdepthrow.ag, genoclonewide)
	  
# Pull out year/season etc
	
	m$clone <- as.character(m$clone)
        temp <- unlist(strsplit(m$clone, split="_"))
        mat <- matrix(temp, ncol=4, byrow=TRUE)
        matdat <- as.data.table(mat)
        m$pop <- matdat$V3
        m$year <- matdat$V2
        m$season <- matdat$V1
        m$seasonB <- ifelse(m$year=="2017" & m$season=="Spring", "April", m$season)
	
	save(m, file="hetandmedrdforSCA_20190506.Rdata")

			
```
