#!/usr/bin/env Rscript

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(phangorn)
  library(doMC)
  registerDoMC(20)
  library(ggplot2)
  library(SNPRelate)
  library(ggbeeswarm)
  #library(tidyverse)


### open genotype file
  genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds",
                      allow.duplicate=TRUE)

### Load LD pruned SNP file to use for Filtering
  filtsnptb <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexandObtusaC/finalsetsnpset01pulex_table_20200207")
  colnames(filtsnptb) <- c("oldvariantids", "chr", "pos", "olddp")

###Optional use non LD pruned dataset
  #filtsnptb <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexandObtusaC/snpsvarpulexpresentinhalf_table_20200207")
  #colnames(filtsnptb) <- c("oldvariantids", "chr", "pos", "olddp")

### Filter SNPs
  snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
      chr = seqGetData(genofile, "chromosome"),
      pos = seqGetData(genofile, "position"))
  setkey(snps, chr, pos)
  setkey(filtsnptb, chr, pos)
  msnps <- merge(filtsnptb, snps)
  msnps[,final.use:=T]

### Load superclone file (does not contain Pulicaria)
  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexandObtusaC/CloneInfoFilePulexandObtusa_withmedrd_update20200324")

### Do some filtering on read depth and independence
  scrd5 <- sc[medrd>4 & Nonindependent==0 & Species=="pulex"]
  #scrd5 <- sc[medrd>9 & Nonindependent==0 & Species=="pulex" & population=="D8"]
  subscrd5 <- scrd5[, c("clone", "SC", "population", "year", "Sex", "Species", "medrd"), with=FALSE]
  clonestouse <- subscrd5
  clonestouseids <- clonestouse$clone

### Optional subsetting of superclones
  #others <- subscrd5[SC=="OO"]
  #scs <- subscrd5[SC!="OO"]
  #listofsc <- unique(scs$SC)
  #setkey(scs, SC)

  #scsub <- foreach(sc.i=listofsc, .combine="rbind")%dopar%{
  #    #sc.i <- "A"
  #  sctmp <- scs[J(sc.i)]
  #  sctmpmrd <- mean(sctmp$medrd)
  #  sctmpB <- sctmp[medrd > sctmpmrd]
  #  sctmpsub <- if(dim(sctmpB)[1] > 1) {sctmpB[sample(nrow(sctmpB), 1),]} else {sctmpB}
  #  sctmpsub
  #  }

  #pulextouse <- rbind(scsub, others)

#  clonestouse <- rbind(subscrd5, pulicaria)
#  clonestouse <- rbind(pulextouse, pulicaria)
  #clonestouse <- pulextouse
  #clonestouseids <- clonestouse$clone

  #seqSetFilter(genofile, variant.id=msnps$variant.ids, sample.id=clonestouse[!is.na(population)]$clone)
  #seqGDS2VCF(genofile, "kinshipmedrd10DBunkSubset.vcf.gz")

### make windows
  setkey(msnps, chr, pos)

  chr.ag <- msnps[,list(start = min(pos), stop=max(pos)), list(chr)]
  chr.ag$length <- chr.ag$stop-chr.ag$start
  chrtouse <- chr.ag$chr[chr.ag$length>1000000]


  #window.size <- 1000000
  #step.size <- 1000000
  window.size <- 1000000
  step.size <- 100000
  wins <- foreach(chr.i=chrtouse, .combine="rbind")%dopar%{
    #chr.i <- "Scaffold_1863_HRSCAF_2081"
    snp.dt.ag <- msnps[J(chr.i)][(final.use), list(start=min(pos), stop=max(pos))]

    data.table(chr=chr.i, start=seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size),
                          stop =seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size) + window.size)
  }

  ### Maybe look at how many SNPs go into each window/tree? Perhaps there are some low ones that are throwing things off?
  msnpcount <- foreach(win.i = c(1:dim(wins)[1]), .combine="rbind", .errorhandling="remove")%dopar%{

    print(paste("Getting data: ", win.i, sep=""))

    snp.dt.tmp <- msnps[J(wins[win.i]$chr)][pos>=wins[win.i]$start & pos<=wins[win.i]$stop][(final.use)]

    tmp <- data.table(window=win.i, numsnps=dim(snp.dt.tmp)[1])
    tmp

  }

  wins$numsnps <- msnpcount$numsnps
  winstouse <- wins[numsnps>500]


### run windows
  m <- foreach(win.i = c(1:dim(winstouse)[1]), .combine="rbind", .errorhandling="remove")%dopar%{
  #m <- foreach(win.i =1:5, .errorhandling="remove", .combine="rbind")%dopar%{

    print(paste("Getting data: ", win.i, sep=""))

    snp.dt.tmp <- msnps[J(winstouse[win.i]$chr)][pos>=winstouse[win.i]$start & pos<=winstouse[win.i]$stop][(final.use)]

    seqResetFilter(genofile)
    seqSetFilter(genofile, variant.id=snp.dt.tmp$variant.ids, sample.id=clonestouse[!is.na(population)]$clone)

    #First lets try IBS with we keep SNPs present in 1 individuals
		### set some global parameters
				maf <- 0.001
				missing.rate <- 0.15
				threads <- 10

				ibs <- snpgdsIBS(genofile, snp.id=snp.dt.tmp$variant.ids, sample.id=clonestouseids, num.thread=20, maf=maf,
							missing.rate=0.15, autosome.only = FALSE)

			### a bit of re-formating of the ibs matrix
				ibs.mat <- ibs$ibs
				rownames(ibs.mat) <- ibs$sample.id
				colnames(ibs.mat) <- ibs$sample.id

			### make the IBs matrix long form
        ibs.matdt <- as.data.table(ibs.mat)
        setkey(clonestouse, clone)
        ibs.matdt$cloneA <- clonestouse$clone
        ibs.long<- melt(ibs.matdt, measure.vars=clonestouseids, variable.name="cloneB", value.name="IBS")
				ibs.long <- na.omit(ibs.long)

        # First let's remove all identical comparisons from ibs.long
  				ibs.longnoident <- ibs.long[ibs.long$cloneA!=ibs.long$cloneB]
          ibs.longnoident$cloneA <- as.factor(ibs.longnoident$cloneA)

  			#Let's also remove duplicated comparisons - how to do this?
  				ibs.longnoident$cloneAnum <- as.numeric(ibs.longnoident$cloneA)
  				ibs.longnoident$cloneBnum <- as.numeric(ibs.longnoident$cloneB)
  				ibs.longnoident$CloneNumComb <- ifelse(ibs.longnoident$cloneAnum > ibs.longnoident$cloneBnum,
  					paste(ibs.longnoident$cloneAnum,ibs.longnoident$cloneBnum,sep="_"),
  					paste(ibs.longnoident$cloneBnum,ibs.longnoident$cloneAnum,sep="_"))
  				setkey(ibs.longnoident, CloneNumComb)
  				ibs.longnoidentunique <- unique(ibs.longnoident, by="CloneNumComb")
  			#Now get back to the original three columns
  				ibs.longunique <- data.table(cloneA=ibs.longnoidentunique$cloneA,
  					cloneB=ibs.longnoidentunique$cloneB, IBS=ibs.longnoidentunique$IBS)

        ibs.longunique$chr <- c(winstouse[win.i]$chr)
        ibs.longunique$start <- c(winstouse[win.i]$start)
        ibs.longunique$stop <- c(winstouse[win.i]$stop)
        ibs.longunique$window <- c(win.i)
        ibs.longunique

  }

  save(m, file="m_IBSbyslidingwindow_subsetSC_1000000_100000_nosub_202003025.Rdata")

  #!/usr/bin/env Rscript

  ### libraries
    library(data.table)
    library(foreach)
    library(ggplot2)

  load("m_IBSbyslidingwindow_subsetSC_1000000_100000_nosub_202003025.Rdata")

  temp <- unlist(strsplit(as.character(m$cloneA), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  m$populationA <- matdat$V3
  m$yearA <- matdat$V2

  temp <- unlist(strsplit(as.character(m$cloneB), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  m$populationB <- matdat$V3
  m$yearB <- matdat$V2

#Now let's add in SC info
  scsubA <- data.table(cloneA=sc$clone, SCA=sc$SC)
  scsubB <- data.table(cloneB=sc$clone, SCB=sc$SC)
  colnames(scsubA) <- c("cloneA", "SCA")
  setkey(scsubA, cloneA)
  setkey(scsubB, cloneB)
  setkey(m, cloneB)
  mtmp <- merge(m, scsubB)
  setkey(mtmp, cloneA)
  msc <- merge(mtmp, scsubA)
  save(msc, file="m_IBSbyslidingwindow_subsetSC_1000000_100000_nosub_202003025_withsc.Rdata")


  msc <- msc[populationA!="DOil" & populationA!="Dramp" & populationB!="DOil" & populationB!="Dramp"]

  msc$populationA <- str_replace(msc$populationA, "Dcat", "DCat")
  msc$populationB <- str_replace(msc$populationB, "Dcat", "DCat")
  msc$pondcompare <- paste(msc$populationA, msc$populationB, sep="_")
  msc$SCcompare <- paste(msc$SCA, msc$SCB, sep="_")

  msc$pondcompare <- str_replace(msc$pondcompare, "D8_D10", "D10_D8")
  msc$pondcompare <- str_replace(msc$pondcompare, "DBunk_D10", "D10_DBunk")
  msc$pondcompare <- str_replace(msc$pondcompare, "DCat_D10", "D10_DCat")
  msc$pondcompare <- str_replace(msc$pondcompare, "DBunk_D8", "D8_DBunk")
  msc$pondcompare <- str_replace(msc$pondcompare, "DCat_D8", "D8_DCat")
  msc$pondcompare <- str_replace(msc$pondcompare, "W1_D8", "D8_W1")
  msc$pondcompare <- str_replace(msc$pondcompare, "W6_D8", "D8_W6")
  msc$pondcompare <- str_replace(msc$pondcompare, "W1_DBunk", "DBunk_W1")
  msc$pondcompare <- str_replace(msc$pondcompare, "W6_DBunk", "DBunk_W6")

  save(msc, file="m_IBSbyslidingwindow_subsetSC_1000000_100000_nosub_202003025_withscsubpond.Rdata")



  msc.ag <- msc[,list(meanIBS=mean(distance), sd_IBS=sd(distance)),
    list(pondcompare, window, chr, start, stop)]
  ABsub <- data.table(pondcompare=c("A_C"), window=msc$window[msc$SCcompare=="A_C"],
    chr=msc$chr[msc$SCcompare=="A_C"], start=msc$start[msc$SCcompare=="A_C"],
    stop=msc$stop[msc$SCcompare=="A_C"], meanIBS=msc$distance[msc$SCcompare=="A_C"],
    sd_IBS=c("NA"))
  msc.ag2 <- rbind(msc.ag, ABsub)
  msc.ag2$category <- ifelse(msc.ag2$pondcompare=="D8_DBunk" | msc.ag2$pondcompare=="D8_DCat" |
    msc.ag2$pondcompare=="DBunk_DCat", "D8surround", ifelse(msc.ag2$pondcompare=="D10_D8" |
    msc.ag2$pondcompare=="D10_DBunk" | msc.ag2$pondcompare=="D10_DCat", "D10vsD8Surround", ifelse(
    msc.ag2$pondcompare=="D8_W1" | msc.ag2$pondcompare=="D8_W6" |
    msc.ag2$pondcompare=="DBunk_W1" | msc.ag2$pondcompare=="DBunk_W6" |
    msc.ag2$pondcompare=="W1_DCat" | msc.ag2$pondcompare=="W6_DCat", "WvsD8Surround", ifelse(
    msc.ag2$pondcompare=="A_C", "zAvsB", ifelse(msc.ag2$pondcompare=="W1_D10" |
    msc.ag2$pondcompare=="W6_D10" | msc.ag2$pondcompare=="W6_W1", 6, "WithinPond")))))
  ggplot(data=mscsubpondB, aes(x=as.factor(pondcompare), y=meanIBS, color=as.factor(withinpond))) + geom_beeswarm()
  ggplot(data=msc.ag2[category!="6"], aes(x=as.factor(window), y=meanIBS, group=as.factor(pondcompare), color=as.factor(category))) + geom_line()
  ggplot(data=msc.ag2[category!="6" & category!="zAvsB"], aes(x=as.factor(window), y=sd_IBS, group=as.factor(pondcompare), color=as.factor(category))) + geom_line()













  setkey(msc, pondcompare, SCcompare, window)
  msc.ag <- msc[,list(meanIBS = mean(distance)),
    list(pondcompare, SCcompare, window, chr, start, stop, populationA, populationB, SCA, SCB)]
  ggplot(data=mscsubpondB, aes(x=window, y=meanIBS, color=as.factor(withinpond))) + geom_point()
  ggplot(data=mscsubpondB[window>0 & window < 11], aes(x=as.factor(withinpond), y=meanIBS,
    color=as.factor(withinpond))) + geom_beeswarm() + facet_wrap(~window)
  mscsubpondB$windowwithin <- paste(mscsubpondB$window, mscsubpondB$withinpond, sep="_")
  ggplot(data=mscsubpondB, aes(x=windowwithin, y=meanIBS, color=as.factor(withinpond))) + geom_violin()
  mscsubpondB$A_C <- ifelse(mscsubpondB$SCcompare=="A_C" | mscsubpondB$SCcompare=="C_A", 1, 0)
  mscsubpondB$withinpondB <- ifelse(mscsubpondB$A_C==1, 2, ifelse(mscsubpondB$A_C==0 & mscsubpondB$withinpond==1, 1, 0))
  ggplot(data=mscsubpondB, aes(x=windowwithin, y=meanIBS, color=as.factor(withinpondB))) + geom_point()
  mscsubpondB$windowwithinB <- paste(mscsubpondB$window, mscsubpondB$withinpondB, sep="_")
  ggplot(data=mscsubpondB, aes(x=windowwithinB, y=meanIBS, color=as.factor(withinpondB))) + geom_point()
  ggplot(data=mscsubpondB[window < 60 & populationA!="D10" & populationB!="D10"], aes(x=windowwithinB, y=meanIBS, color=as.factor(withinpondB))) + geom_point()

  msc.ag$pondcompare <- str_replace(msc.ag$pondcompare, "D8_D10", "D10_D8")
  msc.ag$pondcompare <- str_replace(msc.ag$pondcompare, "DBunk_D10", "D10_DBunk")
  msc.ag$pondcompare <- str_replace(msc.ag$pondcompare, "DCat_D10", "D10_DCat")
  msc.ag$pondcompare <- str_replace(msc.ag$pondcompare, "DBunk_D8", "D8_DBunk")
  msc.ag$pondcompare <- str_replace(msc.ag$pondcompare, "DCat_D8", "D8_DCat")
  msc.ag$pondcompare <- str_replace(msc.ag$pondcompare, "DCat_DBunk", "DBunk_DCat")

  msc.ag2 <- msc.ag[,list(meanIBS=mean(meanIBS), sd_IBS=sd(meanIBS)),
    list(pondcompare, window, chr, start, stop, populationA, populationB)]

  mscsubpond <- msc.ag2[populationA=="D8" | populationA=="DBunk" | populationA=="DCat" | populationA=="D10"]
  mscsubpondB <- mscsubpond[populationB=="D8" | populationB=="DBunk" | populationB=="DCat" | populationB=="D10"]
  mscsubpondB$withinpond <- ifelse(mscsubpondB$populationA==mscsubpondB$populationB, 1, 0)
  ggplot(data=mscsubpondB, aes(x=as.factor(pondcompare), y=meanIBS, color=as.factor(withinpond))) + geom_beeswarm()
