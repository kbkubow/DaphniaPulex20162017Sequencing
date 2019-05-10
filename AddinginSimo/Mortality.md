Let's look at mortality over time, and how this varies for different clones. Overall there are some issues with this.
We only know clonal lineage for clones we sequenced, and they had to survive long enough in the lab and establish well enough to be sequenced.
So aleady, there is bias in who was sequenced. There was a lot of early mortality when the clones arrived in the lab, and we don't know what clonal lineage those clones were.
But let's look at what we do have. This is for the 2017 clones, since we have sequence for ~100 of them for D8 and DBunk.
This is split between individuals that arrived in April, and those at the beginning of May.
```
#!/usr/bin/env Rscript

### libraries
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(tidyverse)

#Read in files
        mort <- fread("D82017ClonesMortality.csv")
        scs <- fread("Superclones20161718withlowcoverageindupdated_20190501")
        
        colnames(mort) <- c("Clone", "DateMort", "Tossed", "AliveMarch7", "DateArr", "DayMort", "Season")
        mort$clone <- paste(mort$Season, "2017", mort$Clone, sep="_")
        mort$clone <- str_replace(mort$clone, "k ", "k_")

        temp <- unlist(strsplit(scs$clone, split="_"))
		    mat <- matrix(temp, ncol=4, byrow=TRUE)
		    matdat <- as.data.table(mat)
	    	scs$population <- matdat$V3
		    scs$year <- matdat$V2
		    scs$season <- matdat$V1
		    scs$clone <- matdat$V4
		    scs$season <- ifelse(scs$year=="2017" & scs$season=="Spring", "April", scs$season)
		    scs$clone <- paste(scs$season, scs$year, scs$population, scs$clone, sep="_")

        setkey(scs, clone)
        setkey(mort, clone)
        m <- merge(scs, mort)
               
        m$alive130day <- ifelse(m$DayMort<130, 0, 1)
        m[is.na(alive130day),alive130day:=1]
        
        m$alive160day <- ifelse(m$DayMort<160, 0, 1)
        m[is.na(alive160day),alive160day:=1]
        
        m$alive190day <- ifelse(m$DayMort<190, 0, 1)
        m[is.na(alive190day),alive190day:=1]
        
        m$alive300day <- ifelse(m$DayMort<300, 0, 1)
        m[is.na(alive300day),alive300day:=1]
        
        m$alive330day <- ifelse(m$DayMort<330, 0, 1)
        m[is.na(alive330day),alive330day:=1]
        
        tmp130 <- as.data.table(table(m$SC, m$alive130))
        colnames(tmp130) <- c("SC", "Alive", "N")
        tmp130wide <- dcast(tmp130, SC ~ Alive, value.var="N")
        colnames(tmp130wide) <- c("SC", "Dead", "Alive")
        tmp130wide$PropAlive <- tmp130wide$Alive/(tmp130wide$Alive+tmp130wide$Dead)
        tmp130wide$Days <- c("130")
 
        tmp160 <- as.data.table(table(m$SC, m$alive160))
        colnames(tmp160) <- c("SC", "Alive", "N")
        tmp160wide <- dcast(tmp160, SC ~ Alive, value.var="N")
        colnames(tmp160wide) <- c("SC", "Dead", "Alive")
        tmp160wide$PropAlive <- tmp160wide$Alive/(tmp160wide$Alive+tmp160wide$Dead)
        tmp160wide$Days <- c("160")

       	tmp190 <- as.data.table(table(m$SC, m$alive190))
        colnames(tmp190) <- c("SC", "Alive", "N")
        tmp190wide <- dcast(tmp190, SC ~ Alive, value.var="N")
        colnames(tmp190wide) <- c("SC", "Dead", "Alive")
        tmp190wide$PropAlive <- tmp190wide$Alive/(tmp190wide$Alive+tmp190wide$Dead)
        tmp190wide$Days <- c("190")

       	tmp300 <- as.data.table(table(m$SC, m$alive300))
        colnames(tmp300) <- c("SC", "Alive", "N")
        tmp300wide <- dcast(tmp300, SC ~ Alive, value.var="N")
        colnames(tmp300wide) <- c("SC", "Dead", "Alive")
        tmp300wide$PropAlive <- tmp300wide$Alive/(tmp300wide$Alive+tmp300wide$Dead)
        tmp300wide$Days <- c("300")

       	tmp330 <- as.data.table(table(m$SC, m$alive330))
        colnames(tmp330) <- c("SC", "Alive", "N")
        tmp330wide <- dcast(tmp330, SC ~ Alive, value.var="N")
        colnames(tmp330wide) <- c("SC", "Dead", "Alive")
        tmp330wide$PropAlive <- tmp330wide$Alive/(tmp330wide$Alive+tmp330wide$Dead)
        tmp330wide$Days <- c("330")

	propmort <- rbind(tmp130wide, tmp160wide, tmp190wide, tmp300wide, tmp330wide)
	propmort$DaysB <- as.integer(propmort$Days)

	ggplot(data=propmort, aes(x=DaysB, y=PropAlive, group=SC, color=SC)) + geom_line()
	ggplot(data=propmort[propmort$Dead+propmort$Alive > 4], aes(x=DaysB, y=PropAlive, group=SC, color=SC)) + geom_line()
	
	scssub <- data.table(SC=scs$SC, population=scs$population)
	scssubuniq <- unique(scssub)
	scssubuniqnoOO <- scssubuniq[SC!="OO"]
	propmortnoOO <- propmort[SC!="OO"]
	setkey(scssubuniqnoOO, SC)
	setkey(propmortnoOO, SC)
	mpropmort <- merge(propmortnoOO, scssubuniqnoOO)
	mpropmort$notuse <- ifelse(mpropmort$SC=="A" & mpropmort$population=="DBunk", 0, ifelse(
		mpropmort$SC=="B" & mpropmort$population=="Dmud", 0, ifelse(mpropmort$SC=="AD", 0, 1)))
	ggplot(data=mpropmort, aes(x=DaysB, y=PropAlive, group=SC, color=SC)) + geom_line()
	ggplot(data=mpropmort[Dead+Alive > 4 & notuse==1], aes(x=DaysB, y=PropAlive, group=SC, color=SC)) + geom_line() + facet_wrap(~population)
	ggplot(data=mpropmort[notuse==1], aes(x=DaysB, y=PropAlive, group=SC, color=SC)) + geom_line() + facet_wrap(~population)
	
	library(lme4)
	mpropmort$total <- mpropmort$Dead+mpropmort$Alive
	mpropmort5 <- mpropmort[total > 4]
	mpropmort5D8 <- mpropmort5[population=="D8"]
	mpropmort5DBunk <- mpropmort5[population=="DBunk" & notuse==1]
	test <- glm(PropAlive~SC+Days+SC*Days, data=mpropmort5, family="binomial", weights=mpropmort5$total)
	test <- glm(PropAlive~SC+DaysB+SC*DaysB, data=mpropmort5D8, family="binomial", weights=mpropmort5D8$total)
	test <- glm(PropAlive~SC+Days+SC*Days, data=mpropmort5D8, family="binomial", weights=mpropmort5D8$total)
	test <- glm(PropAlive~SC*DaysB, data=mpropmort5D8, family="binomial", weights=mpropmort5D8$total)

	mpropmort5D8330 <- mpropmort5D8[Days=="330"]
	test330 <- glm(PropAlive~SC, data=mpropmort5D8330, family="binomial", weights=mpropmort5D8330$total)
	summary(test330)
	mpropmort5D8300 <- mpropmort5D8[Days=="300"]
	test300 <- glm(PropAlive~SC, data=mpropmort5D8300, family="binomial", weights=mpropmort5D8300$total)
	summary(test300)
	mpropmort5D8190 <- mpropmort5D8[Days=="190"]
	test190 <- glm(PropAlive~SC, data=mpropmort5D8190, family="binomial", weights=mpropmort5D8190$total)
	summary(test190)
	mpropmort5D8160 <- mpropmort5D8[Days=="160"]
	test160 <- glm(PropAlive~SC, data=mpropmort5D8160, family="binomial", weights=mpropmort5D8160$total)
	summary(test160)
	mpropmort5D8130 <- mpropmort5D8[Days=="130"]
	test130 <- glm(PropAlive~SC, data=mpropmort5D8130, family="binomial", weights=mpropmort5D8130$total)
	summary(test130)


	mpropmort5DBunk330 <- mpropmort5DBunk[Days=="330"]
	test330Bunk <- glm(PropAlive~SC, data=mpropmort5DBunk330, family="binomial", weights=mpropmort5DBunk330$total)
	summary(test330Bunk)
	mpropmort5DBunk300 <- mpropmort5DBunk[Days=="300"]
	test300Bunk <- glm(PropAlive~SC, data=mpropmort5DBunk300, family="binomial", weights=mpropmort5DBunk300$total)
	summary(test300Bunk)
	mpropmort5DBunk190 <- mpropmort5DBunk[Days=="190"]
	test190Bunk <- glm(PropAlive~SC, data=mpropmort5DBunk190, family="binomial", weights=mpropmort5DBunk190$total)
	summary(test190Bunk)


	m$DayMortB <- m$DayMort
	m[is.na(DayMortB),DayMortB:=350]
	mtouse <- m[SC=="A" | SC=="B" | SC=="K" | SC=="D" | SC=="E" | SC=="F" | SC=="J" | SC=="K"]
	survdiff(DayMortB~SC, mtouse)
