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

# Pull out all A, B, M, amd DCat18004 clones

        scsABM <- scs[superClone.index=="2" | superClone.index=="3" | superClone.index=="33" | superClone.index=="86"]

#Remove males and SM libraries

        scsAnomaleSM <- scsABM[clone!="Lab_2019_D8_349Male" & clone!="May_2017_D8_770SM"]
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
              
              setkey(mhetlong, clone)
              setkey(scs, clone)
              mmhetlong <- merge(mhetlong, scs)
              
              mmhetlongA <- mmhetlong[SC=="A"]            

#Count instances of variant.ids/clone/dosage	
              genoclone <- mmhetlongA[, .N, by=list(variant.ids, dosage)]

#Remove NAs
	      genoclone <- genoclone[dosage!="NA"]
              
#Transform to wide format
	genoclonewide <- dcast(genoclone, variant.ids ~ dosage, value.var="N")
	colnames(genoclonewide) <- c("variant.id", "dos0", "dos1", "dos2")
        genoclonewide[is.na(dos0),dos0:=0]
        genoclonewide[is.na(dos1),dos1:=0]
        genoclonewide[is.na(dos2),dos2:=0]

        genoclonewide$consensusdosage <- ifelse(genoclonewide$dos2 > genoclonewide$dos1 & 
                genoclonewide$dos2 > genoclonewide$dos0, "2", ifelse(genoclonewide$dos0 > genoclonewide$dos1 & 
                genoclonewide$dos0 > genoclonewide$dos2, "0", "1"))

        genoclonewidesubA <- data.table(variant.id=genoclonewide$variant.id, Acondosage=genoclonewide$consensusdosage)

              mmhetlongB <- mmhetlong[SC=="B"]            

#Count instances of variant.ids/clone/dosage	
              genoclone <- mmhetlongB[, .N, by=list(variant.ids, dosage)]

#Remove NAs
	      genoclone <- genoclone[dosage!="NA"]
              
#Transform to wide format
	genoclonewide <- dcast(genoclone, variant.ids ~ dosage, value.var="N")
	colnames(genoclonewide) <- c("variant.id", "dos0", "dos1", "dos2")
        genoclonewide[is.na(dos0),dos0:=0]
        genoclonewide[is.na(dos1),dos1:=0]
        genoclonewide[is.na(dos2),dos2:=0]

        genoclonewide$consensusdosage <- ifelse(genoclonewide$dos2 > genoclonewide$dos1 & 
                genoclonewide$dos2 > genoclonewide$dos0, "2", ifelse(genoclonewide$dos0 > genoclonewide$dos1 & 
                genoclonewide$dos0 > genoclonewide$dos2, "0", "1"))

        genoclonewidesubB <- data.table(variant.id=genoclonewide$variant.id, Acondosage=genoclonewide$consensusdosage)
        colnames(genoclonewidesubB) <- c("variant.id", "Bcondosage")

              mmhetlongM <- mmhetlong[SC=="M"]            

#Count instances of variant.ids/clone/dosage	
              genoclone <- mmhetlongM[, .N, by=list(variant.ids, dosage)]

#Remove NAs
	      genoclone <- genoclone[dosage!="NA"]
              
#Transform to wide format
	genoclonewide <- dcast(genoclone, variant.ids ~ dosage, value.var="N")
	colnames(genoclonewide) <- c("variant.id", "dos0", "dos1", "dos2")
        genoclonewide[is.na(dos0),dos0:=0]
        genoclonewide[is.na(dos1),dos1:=0]
        genoclonewide[is.na(dos2),dos2:=0]

        genoclonewide$consensusdosage <- ifelse(genoclonewide$dos2 > genoclonewide$dos1 & 
                genoclonewide$dos2 > genoclonewide$dos0, "2", ifelse(genoclonewide$dos0 > genoclonewide$dos1 & 
                genoclonewide$dos0 > genoclonewide$dos2, "0", "1"))

        genoclonewidesubM <- data.table(variant.id=genoclonewide$variant.id, Acondosage=genoclonewide$consensusdosage)
        colnames(genoclonewidesubM) <- c("variant.id", "Mcondosage")

              mmhetlongO <- mmhetlong[SC=="OO"]            

        genoclonewidesubO <- data.table(variant.id=mmhetlongO$variant.ids, Ocondosage=mmhetlongO$dosage)
        
        setkey(genoclonewidesubA, variant.id)
        setkey(genoclonewidesubB, variant.id)
        setkey(genoclonewidesubM, variant.id)
        setkey(genoclonewidesubO, variant.id)
        mcon <- merge(genoclonewidesubA, genoclonewidesubB)
        setkey(mcon, variant.id)
        mmcon <- merge(mcon, genoclonewidesubM)
        setkey(mmcon, variant,id)
        mmmcon <- merge(mmcon, genoclonewidesubO)
       
	save(mmmcon, file="mmmcon_20190723.Rdata")


	scst <- data.table(sc=c("Acondosage", "Bcondosage", "Mcondosage", "Ocondosage"))
		
	scsfull <- foreach(x=(1:4), .combine="rbind")%do%{
		
		SCA <- scst$sc[x]
		data.table(SCA=c(SCA), SCB=scst$sc)
		
		}
	
	scsfullnoident <- scsfull[SCA!=SCB]
	
	save(scsfullnoident, file="scsfullnoident_20190723.Rdata")




```
