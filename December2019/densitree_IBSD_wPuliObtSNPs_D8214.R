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
genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")

###Use non LD pruned dataset
  load("dpfiltsnps_20200224.Rdata")
  filtsnptb <- dpfiltsnps
  colnames(filtsnptb) <- c("oldvariantids", "chr", "pos", "olddp")

### Filter SNPs
  snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
      chr = seqGetData(genofile, "chromosome"),
      pos = seqGetData(genofile, "position"))
  setkey(snps, chr, pos)
  setkey(filtsnptb, chr, pos)
  msnps <- merge(filtsnptb, snps)
  msnps[,final.use:=T]

### Load superclone file
  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/Superclones201617182019withObtusaandPulicaria_kingcorr_20200402_wmedrd.txt")

### Do some filtering on read depth and independence
  scrd5 <- sc[medrd>4 & Nonindependent==0]
  subscrd5 <- scrd5[, c("clone", "SC", "population", "year", "Sex", "Species", "medrd"), with=FALSE]
  clonestouse <- subscrd5
  clonestouseids <- clonestouse$clone

### make windows
  setkey(msnps, chr, pos)

  chr.ag <- msnps[,list(start = min(pos), stop=max(pos)), list(chr)]
  chr.ag$length <- chr.ag$stop-chr.ag$start
  chrtouse <- chr.ag$chr[chr.ag$length>1000000]

  window.size <- 250000
  step.size <- 10000
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
  winstouse <- wins[numsnps>1000]


### run windows
  m <- foreach(win.i = c(1:dim(winstouse)[1]), .combine="rbind", .errorhandling="remove")%dopar%{
  #m <- foreach(win.i =varA:varB, .errorhandling="remove", .combine="rbind")%dopar%{

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
          tmpnum <- data.table(cloneB=ibs.longnoident$cloneA, cloneBnum=ibs.longnoident$cloneAnum)
          tmpnumu <- unique(tmpnum)
          setkey(ibs.longnoident, cloneB)
          setkey(tmpnumu, cloneB)
          mibs.longnoident <- merge(ibs.longnoident, tmpnumu)
  				mibs.longnoident$CloneNumComb <- ifelse(mibs.longnoident$cloneAnum > mibs.longnoident$cloneBnum,
  					paste(mibs.longnoident$cloneAnum,mibs.longnoident$cloneBnum,sep="_"),
  					paste(mibs.longnoident$cloneBnum,mibs.longnoident$cloneAnum,sep="_"))
  				setkey(mibs.longnoident, CloneNumComb)
  				ibs.longnoidentunique <- unique(mibs.longnoident, by="CloneNumComb")
  			#Now get back to the original three columns
  				ibs.longunique <- data.table(cloneA=ibs.longnoidentunique$cloneA,
  					cloneB=ibs.longnoidentunique$cloneB, IBS=ibs.longnoidentunique$IBS)

        ibs.longunique$chr <- c(winstouse[win.i]$chr)
        ibs.longunique$start <- c(winstouse[win.i]$start)
        ibs.longunique$stop <- c(winstouse[win.i]$stop)
        ibs.longunique$window <- c(win.i)
        ibs.longunique

  }

  save(m, file="m_IBSbyslidingwindow_250000_10000_withpulicariaandobtusa_202003030.Rdata")

  #!/usr/bin/env Rscript

  ### libraries
    library(data.table)
    library(foreach)
    library(ggplot2)

  ### Let's work through one input file first
    load(paste("m_IBSbyslidingwindow_250000_10000_withpulicariaandobtusa_20200424_", varA, ".Rdata", sep=""))

  # Load superclone file
    sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/Superclones201617182019withObtusaandPulicaria_kingcorr_20200402_wmedrd.txt")

  # Add in superclone information
    scsubA <- data.table(cloneA=sc$clone, SCA=sc$SC, popA=sc$population, yearA=sc$year, medrdA=sc$medrd, speciesA=sc$Species, sexA=sc$Sex)
    scsubB <- data.table(cloneB=sc$clone, SCB=sc$SC, popB=sc$population, yearB=sc$year, medrdB=sc$medrd, speciesB=sc$Species, sexB=sc$Sex)
    setkey(scsubA, cloneA)
    setkey(scsubB, cloneB)
    setkey(m, cloneB)
    mtmp <- merge(m, scsubB)
    setkey(mtmp, cloneA)
    msc <- merge(mtmp, scsubA)


  # Pull out just D8 and DBunk
    mscD8DBunk <- msc[popA=="D8" | popA=="DBunk"]
    mscD8DBunkB <- mscD8DBunk[popB=="D8" | popB=="DBunk"]
    mscD8DBunkB$SCcomp <- paste(mscD8DBunkB$SCA, mscD8DBunkB$SCB, sep="_")
    mscD8DBunkBAC <- mscD8DBunkB[SCcomp=="A_C" | SCcomp=="C_A" | SCcomp=="A_A" | SCcomp=="C_C"]

    mscD8DBunkBACD8214vsC <- mscD8DBunkBAC[cloneA=="April_2017_D8_214" & SCB=="C" | cloneB=="April_2017_D8_214" & SCA=="C"]
    mscD8DBunkBACD8214vsA <- mscD8DBunkBAC[cloneA=="April_2017_D8_214" & SCB=="A" | cloneB=="April_2017_D8_214" & SCA=="A"]
    mscD8DBunkBACD8823vsC <- mscD8DBunkBAC[cloneA=="Spring_2016_D8_8.23" & SCB=="C" | cloneB=="Spring_2016_D8_8.23" & SCA=="C"]
    mscD8DBunkBACD8823vsA <- mscD8DBunkBAC[cloneA=="Spring_2016_D8_8.23" & SCB=="A" | cloneB=="Spring_2016_D8_8.23" & SCA=="A"]
    mscD8DBunkBACD8823vsD8214 <- mscD8DBunkBAC[cloneA=="Spring_2016_D8_8.23" & cloneB=="April_2017_D8_214" |
      cloneB=="Spring_2016_D8_8.23" & cloneA=="April_2017_D8_214"]

  # Now aggregate and combine
    mscD8DBunkBACD8214vsC.ag <- mscD8DBunkBACD8214vsC[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
      list(window, chr, start, stop)]
      mscD8DBunkBACD8214vsC.ag$compinfo <- c("D8214vsC")
    mscD8DBunkBACD8214vsA.ag <- mscD8DBunkBACD8214vsA[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
      list(window, chr, start, stop)]
      mscD8DBunkBACD8214vsA.ag$compinfo <- c("D8214vsA")
    mscD8DBunkBACD8823vsC.ag <- mscD8DBunkBACD8823vsC[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
      list(window, chr, start, stop)]
      mscD8DBunkBACD8823vsC.ag$compinfo <- c("D8823vsC")
    mscD8DBunkBACD8823vsA.ag <- mscD8DBunkBACD8823vsA[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
      list(window, chr, start, stop)]
      mscD8DBunkBACD8823vsA.ag$compinfo <- c("D8823vsA")
    mscD8DBunkBACD8823vsD8214.ag <- mscD8DBunkBACD8823vsD8214[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
      list(window, chr, start, stop)]
      mscD8DBunkBACD8823vsD8214.ag$compinfo <- c("D8823vsD8214")

    mscD8DBunkBACD8823vsD8214tot <- rbind(mscD8DBunkBACD8214vsC.ag, mscD8DBunkBACD8214vsA.ag,
      mscD8DBunkBACD8823vsC.ag, mscD8DBunkBACD8823vsA.ag, mscD8DBunkBACD8823vsD8214.ag)

    save(mscD8DBunkBACD8823vsD8214tot.ag, file=paste("mscD8DBunkBACD8823vsD8214tot", varA, ".Rdata", sep=""))


#!/usr/bin/env Rscript

    ### libraries
      library(data.table)
      library(foreach)
      library(ggplot2)

  ### Load processed IBS by window files
      inputfiles <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", pattern="mscACvsrest.ag_20200424")

      totalibsag <- foreach(i=1:length(inputfiles), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", inputfiles[i], sep="")
        load(f)
        mscACvsrest.ag
      }

      inputfilesPO <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", pattern="mscACvsrest_POBC_20200424")

      totalibsagPO <- foreach(i=1:length(inputfilesPO), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", inputfilesPO[i], sep="")
        load(f)
        mscACvsrest_POBC
      }

      inputfilesDBunkD8 <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", pattern="mscD8DBunknoAC_20200424")

      totalibsagD8DBunk <- foreach(i=1:length(inputfilesDBunkD8), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", inputfilesDBunkD8[i], sep="")
        load(f)
        mscD8DBunkBnoAC.ag
      }


      #inputfilesD8214D823 <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", pattern="mscD8DBunkBACD8823vsD8214tot")

      #totalibsD8214D823 <- foreach(i=1:length(inputfilesD8214D823), .combine="rbind")%do%{
      #  f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", inputfilesD8214D823[i], sep="")
      #  load(f)
      #  mscD8DBunkBACD8823vsD8214tot
      #}

      inputfilesObtusaPulicaria <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", pattern="mscothervsObtusaPulicaria")


      totalibsObPul <- foreach(i=1:length(inputfilesObtusaPulicaria), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/WithPulicaria/", inputfilesObtusaPulicaria[i], sep="")
        load(f)
        mscothervsObtusaPulicaria.ag
      }

      totalibsObPulObonly <- totalibsObPul[speciesA=="obtusa" | speciesB=="obtusa"]
      totalibsObPulObonlyD8DBunk <- totalibsObPulObonly[popA=="D8" | popB=="D8" | popA=="DBunk" | popB=="DBunk"]

      totalibsObPulObonlyD8DBunk$outlier <- ifelse(totalibsObPulObonlyD8DBunk$meanIBS > 0.6, 1, 0)
      totalibsObPulObonlyD8DBunknoDbarborpulicaria <- totalibsObPulObonlyD8DBunk[popA!="Dbarb" & popB!="Dbarb" & speciesA!="pulicaria" & speciesB!="pulicaria"]


      ggplot(data=totalibsObPulObonlyD8DBunknoDbarborpulicaria, aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point() +
        facet_wrap(~chr, scales="free_x")

        ggplot(data=totalibsObPulObonlyD8DBunknoDbarborpulicaria[chr=="Scaffold_6786_HRSCAF_7541"], aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point()
        ggplot(data=totalibsObPulObonlyD8DBunknoDbarborpulicaria[chr=="Scaffold_6786_HRSCAF_7541" & window > 5250 & window < 5430], aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point()
        ggplot(data=totalibsObPulObonlyD8DBunknoDbarborpulicaria[chr=="Scaffold_1863_HRSCAF_2081"], aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point()

      totalibsObPulObonlyDCat <- totalibsObPulObonly[popA=="DCat" | popB=="DCat"]
      totalibsObPulObonlyDCat$outlier <- ifelse(totalibsObPulObonlyDCat$meanIBS > 0.6, 1, 0)
      totalibsObPulObonlyDCatnoDbarborpulicaria <- totalibsObPulObonlyDCat[popA!="Dbarb" & popB!="Dbarb" & speciesA!="pulicaria" & speciesB!="pulicaria"]

      ggplot(data=totalibsObPulObonlyDCatnoDbarborpulicaria, aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point() +
        facet_wrap(~chr, scales="free_x")
      ggplot(data=totalibsObPulObonlyDCatnoDbarborpulicaria[chr=="Scaffold_6786_HRSCAF_7541" & window > 5250 & window < 5430], aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point()

      totalibsObPulObonlyDD10 <- totalibsObPulObonly[popA=="D10" | popB=="D10"]
      totalibsObPulObonlyDD10$outlier <- ifelse(totalibsObPulObonlyDD10$meanIBS > 0.6, 1, 0)
      totalibsObPulObonlyDD10noDbarborpulicaria <- totalibsObPulObonlyDD10[popA!="Dbarb" & popB!="Dbarb" & speciesA!="pulicaria" & speciesB!="pulicaria"]

      ggplot(data=totalibsObPulObonlyDD10noDbarborpulicaria, aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point() +
        facet_wrap(~chr, scales="free_x")
      ggplot(data=totalibsObPulObonlyDD10noDbarborpulicaria[chr=="Scaffold_6786_HRSCAF_7541" & window > 5250 & window < 5430], aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point()

      totalibsObPulObonlyW1 <- totalibsObPulObonly[popA=="W1" | popB=="W1"]
      totalibsObPulObonlyW1$outlier <- ifelse(totalibsObPulObonlyW1$meanIBS > 0.6, 1, 0)
      totalibsObPulObonlyW1noDbarborpulicaria <- totalibsObPulObonlyW1[popA!="Dbarb" & popB!="Dbarb" & speciesA!="pulicaria" & speciesB!="pulicaria"]

      ggplot(data=totalibsObPulObonlyW1noDbarborpulicaria, aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point() +
        facet_wrap(~chr, scales="free_x")
      ggplot(data=totalibsObPulObonlyW1noDbarborpulicaria[chr=="Scaffold_6786_HRSCAF_7541" & window > 5250 & window < 5430], aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point()


      totalibsObPulPulonly <- totalibsObPul[speciesA=="pulicaria" | speciesB=="pulicaria"]
      totalibsObPulPulonlyD8DBunk <- totalibsObPulPulonly[popA=="D8" | popB=="D8" | popA=="DBunk" | popB=="DBunk"]

      totalibsObPulPulonlyD8DBunk$outlier <- ifelse(totalibsObPulPulonlyD8DBunk$meanIBS > 0.6, 1, 0)
      totalibsObPulPulonlyD8DBunknoDbarborobtusa <- totalibsObPulPulonlyD8DBunk[popA!="Dbarb" & popB!="Dbarb" & speciesA!="obtusa" & speciesB!="obtusa"]


      ggplot(data=totalibsObPulPulonlyD8DBunknoDbarborobtusa, aes(x=as.factor(window), y=meanIBS)) + geom_point() +
        facet_wrap(~chr, scales="free_x")

        ggplot(data=totalibsObPulPulonlyD8DBunknoDbarborobtusa[chr=="Scaffold_6786_HRSCAF_7541"], aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point()
        ggplot(data=totalibsObPulPulonlyD8DBunknoDbarborobtusa[chr=="Scaffold_6786_HRSCAF_7541" & window > 5250 & window < 5430], aes(x=as.factor(window), y=meanIBS, color=as.factor(outlier))) + geom_point()


### trial with subset of data
    #load("mscACvsrest.ag_20200424_G.Rdata")
    #load("mscACvsrest_POBC_20200424G.Rdata")
    #load("mscD8DBunknoAC_20200424G.Rdata")

    mscD8DBunkBnoAC.ag <- totalibsagD8DBunk

    colnames(mscD8DBunkBnoAC.ag) <- c("window", "chr", "start", "stop", "SCA_C",
      "SCB_C", "yearA", "yearB", "popA", "popB", "compinfo", "speciesA",
      "speciesB", "meanIBS", "sd_IBS")

    mscD8DBunkBnoAC.ag$ACvspond <- paste(mscD8DBunkBnoAC.ag$popA, mscD8DBunkBnoAC.ag$popB, sep="_")

    mscD8DBunkBnoAC.ag$SCcompare <- mscD8DBunkBnoAC.ag$ACvspond
    mscD8DBunkBnoAC.ag <- mscD8DBunkBnoAC.ag[speciesA!="obtusa" & speciesB!="obtusa"]

    mscD8DBunkBnoAC.agsub <- mscD8DBunkBnoAC.ag[, c("window", "chr", "start", "stop", "SCA_C",
      "SCB_C", "ACvspond", "meanIBS", "sd_IBS", "speciesA", "speciesB"), with=FALSE]

    mscD8DBunkBnoAC.agsub <- mscD8DBunkBnoAC.agsub[SCA_C!="G" & SCB_C!="G"]

    #mscACvsrest_AW.ag <- rbind(mscACvsrest.ag, mscACvsrest_POBC, mscD8DBunkBnoAC.agsub)
    mscACvsrest_AW.ag <- rbind(totalibsag, totalibsagPO, mscD8DBunkBnoAC.agsub)
    #mscACvsrest_AW.ag <- rbind(totalibsag, totalibsagPO)


    mscACvsrest_AW.ag$SCcompare <- paste(mscACvsrest_AW.ag$SCA_C,mscACvsrest_AW.ag$SCB_C,sep="_")

    mscACvsrest_AW.ag$ACvspond <- ifelse(mscACvsrest_AW.ag$SCA_C=="C" & mscACvsrest_AW.ag$SCB_C=="G", "C_Obtusa", ifelse(
      mscACvsrest_AW.ag$SCA_C=="C" & mscACvsrest_AW.ag$SCB_C=="P", "C_Pulicaria", mscACvsrest_AW.ag$ACvspond))

    mscACvsrest_AW.ag$ACvspond <- ifelse(mscACvsrest_AW.ag$ACvspond=="DBunk_D8", "D8_DBunk", mscACvsrest_AW.ag$ACvspond)

    mscACvsrest_AW.ag$ACvspond <- factor(mscACvsrest_AW.ag$ACvspond, levels=c("A_D8DBunk_C_D8DBunk",
      "D8_D8", "DBunk_DBunk", "D8_DBunk",
      "A_D10", "A_W1", "A_W6", "A_pulicaria", "A_obtusa", "C_D10", "C_W1", "C_W6", "C_Pulicaria", "C_Obtusa"))

    #mscACvsrest_AW.ag$SCcompare <- factor(mscACvsrest_AW.ag$SCcompare, levels=c("A_C", "A_D", "A_AF", "A_OO_Fall_2016_D10_41",
    #  "C_D", "C_AF", "C_OO_Fall_2016_D10_41", "A_OO_Spring_2016_W1_1.1", "A_OO_Spring_2016_W6_6.1", "A_OO_Spring_2016_W6_6.3",
    #  "C_OO_Spring_2016_W1_1.1", "C_OO_Spring_2016_W6_6.1", "C_OO_Spring_2016_W6_6.3", "A_P", "C_P", "A_G", "C_G"))

    mscACvsrest_AW.ag$colortype <- ifelse(mscACvsrest_AW.ag$ACvspond=="A_D8DBunk_C_D8DBunk", "AvsC",
      ifelse(mscACvsrest_AW.ag$ACvspond=="D8_D8" | mscACvsrest_AW.ag$ACvspond=="D8_DBunk" |
      mscACvsrest_AW.ag$ACvspond=="DBunk_DBunk", "D8_DBunk", ifelse(mscACvsrest_AW.ag$ACvspond=="A_D10" |
      mscACvsrest_AW.ag$ACvspond=="A_W1" | mscACvsrest_AW.ag$ACvspond=="A_W6" |
      mscACvsrest_AW.ag$ACvspond=="A_pulicaria" | mscACvsrest_AW.ag$ACvspond=="A_obtusa", "AvsOther",
      ifelse(mscACvsrest_AW.ag$ACvspond=="C_D10" | mscACvsrest_AW.ag$ACvspond=="C_W1" |
      mscACvsrest_AW.ag$ACvspond=="C_W6" | mscACvsrest_AW.ag$ACvspond=="C_Pulicaria" |
      mscACvsrest_AW.ag$ACvspond=="C_Obtusa", "CvsOther", "NA"))))

    mscACvsrest_AW.ag$colortype <- factor(mscACvsrest_AW.ag$colortype, levels=c("AvsC",
      "D8_DBunk", "AvsOther", "CvsOther"))

    #ggplot(data=mscACvsrest_AW.ag, aes(x=SCcompare, y=meanIBS, group=SCcompare, color=colortype)) + geom_boxplot()

    allbox <- ggplot(data=mscACvsrest_AW.ag, aes(x=ACvspond, y=meanIBS, group=ACvspond, color=colortype)) + geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    allboxnoD8DBunk <- ggplot(data=mscACvsrest_AW.ag[colortype!="D8_DBunk"],
      aes(x=ACvspond, y=meanIBS, group=ACvspond, color=colortype)) + geom_boxplot() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggsave(allbox, file="IBSboxwhiskerall_20200524.pdf")
    ggsave(allboxnoD8DBunk, file="IBSboxwhiskernoD8DBunk_20200524.pdf")


    ggplot(data=mscACvsrest_AW.ag[ACvspond=="A_pulicaria" | ACvspond=="A_obtusa" |
      ACvspond=="C_Pulicaria" | ACvspond=="C_Obtusa" | ACvspond=="A_D8DBunk_C_D8DBunk"],
      aes(x=as.factor(window), y=meanIBS, group=ACvspond, color=ACvspond)) + geom_line() + facet_wrap(~chr, scales="free_x")

    tmp <- mscACvsrest_AW.ag[ACvspond=="A_pulicaria" | ACvspond=="A_obtusa" |
      ACvspond=="C_Pulicaria" | ACvspond=="C_Obtusa" | ACvspond=="A_D8DBunk_C_D8DBunk"]

      ggplot(data=tmp[chr=="Scaffold_6786_HRSCAF_7541" & window > 5250 & window < 5430],
        aes(x=as.factor(window), y=meanIBS, color=ACvspond)) + geom_point() + facet_wrap(~chr, scales="free_x")


    mscACvsrest_AW.ag <- mscACvsrest_AW.ag[window > 10000]

    ggplot() + geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS,
        group=SCB_C, color="AvsW6"), alpha=0.5) +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
          group=SCB_C, color="AvsW6"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
          group=SCB_C, color="AvsW6"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS,
          group=SCB_C, color="AvsW6"), alpha=0.5) +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
          group=SCB_C, color="AvsW6"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
          group=SCB_C, color="AvsW6"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS,
          group=SCB_C, color="AvsW1"), alpha=0.5) +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
          group=SCB_C, color="AvsW1"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
          group=SCB_C, color="AvsW1"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS,
          group=SCB_C, color="AvsD10"), alpha=0.5) +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
          group=SCB_C, color="AvsD10"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
          group=SCB_C, color="AvsD10"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS,
          group=SCB_C, color="AvsD10"), alpha=0.5) +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
          group=SCB_C, color="AvsD10"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
          group=SCB_C, color="AvsD10"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS,
          group=SCB_C, color="AvsD10"), alpha=0.5) +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
          group=SCB_C, color="AvsD10"), alpha=0.5, linetype="dotted") +
      geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
          group=SCB_C, color="AvsD10"), alpha=0.5, linetype="dotted") +
            geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS,
                group=SCB_C, color="CvsW6"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="CvsW6"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="CvsW6"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS,
                  group=SCB_C, color="CvsW6"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="CvsW6"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="CvsW6"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS,
                  group=SCB_C, color="CvsW1"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="CvsW1"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="CvsW1"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS,
                  group=SCB_C, color="CvsD10"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="CvsD10"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="CvsD10"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS,
                  group=SCB_C, color="CvsD10"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="CvsD10"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="CvsD10"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS,
                  group=SCB_C, color="CvsD10"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="CvsD10"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="CvsD10"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS,
                  group=SCB_C, color="CvsPulicaria"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="CvsPulicaria"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="CvsPulicaria"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS,
                  group=SCB_C, color="CvsObtusa"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="CvsObtusa"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="CvsObtusa"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS,
                  group=SCB_C, color="AvsObtusa"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="AvsObtusa"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="AvsObtusa"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS,
                  group=SCB_C, color="AvsPulicaria"), alpha=0.5) +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                  group=SCB_C, color="AvsPulicaria"), alpha=0.5, linetype="dotted") +
              geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                  group=SCB_C, color="AvsPulicaria"), alpha=0.5, linetype="dotted") +
                  geom_line(data=mscACvsrest_AW.ag[SCB_C=="C"], aes(x=as.factor(window), y=meanIBS, group=SCB_C, color="AvsC"), size=1) +
                      geom_line(data=mscACvsrest_AW.ag[SCB_C=="C"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS, group=SCB_C,
                        color="AvsC"), linetype="dotted") +
                      geom_line(data=mscACvsrest_AW.ag[SCB_C=="C"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS, group=SCB_C,
                        color="AvsC"), linetype="dotted")


          ACvsOB <-              ggplot() + geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS,
                                      group=SCB_C, color="CvsPulicaria"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="CvsPulicaria"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="CvsPulicaria"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS,
                                      group=SCB_C, color="CvsObtusa"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="CvsObtusa"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="CvsObtusa"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS,
                                      group=SCB_C, color="AvsObtusa"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="AvsObtusa"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="G"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="AvsObtusa"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS,
                                      group=SCB_C, color="AvsPulicaria"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="AvsPulicaria"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="AvsPulicaria"), alpha=0.5, linetype="dotted") +
                                      geom_line(data=mscACvsrest_AW.ag[SCB_C=="C"], aes(x=as.factor(window), y=meanIBS, group=SCB_C, color="AvsC"), size=1) +
                                          geom_line(data=mscACvsrest_AW.ag[SCB_C=="C"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS, group=SCB_C,
                                            color="AvsC"), linetype="dotted") +
                                          geom_line(data=mscACvsrest_AW.ag[SCB_C=="C"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS, group=SCB_C,
                                            color="AvsC"), linetype="dotted") + facet_wrap(~mscACvsrest_AW.ag$chr)


                      ggplot() + geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS,
                              group=SCB_C, color="CvsPulicaria"), alpha=0.5) +
                      geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                              group=SCB_C, color="CvsPulicaria"), alpha=0.5, linetype="dotted") +
                      geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="P"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                              group=SCB_C, color="CvsPulicaria"), alpha=0.5, linetype="dotted") + facet_wrap(~chr, scales="free_x")

      mscACvsrest_AW.agPO <- mscACvsrest_AW.ag[SCB_C=="P" | SCB_C=="G" | SCB_C=="C"]
      mscACvsrest_AW.agPO$ACvspond <- ifelse(mscACvsrest_AW.agPO$SCA_C=="C" & mscACvsrest_AW.agPO$SCB_C=="G", "C_Obtusa", ifelse(
        mscACvsrest_AW.agPO$SCA_C=="C" & mscACvsrest_AW.agPO$SCB_C=="P", "C_Pulicaria", ifelse(mscACvsrest_AW.agPO$SCA_C=="A" &
        mscACvsrest_AW.agPO$SCB_C=="G", "A_Obtusa", ifelse(mscACvsrest_AW.agPO$SCA_C=="A" &
        mscACvsrest_AW.agPO$SCB_C=="P", "A_Pulicaria", mscACvsrest_AW.agPO$ACvspond))))


          ggplot(data=mscACvsrest_AW.agPO, aes(x=as.factor(window), y=meanIBS, group=ACvspond, color=ACvspond)) + geom_line() + facet_wrap(~chr, scales="free_x")
          ggplot(data=mscACvsrest_AW.ag, aes(x=as.factor(window), y=meanIBS, group=SCcompare, color=ACvspond, group=ACvspond)) + geom_line() + facet_wrap(~chr, scales="free_x")

                        ggplot() + geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS,
                            group=SCB_C, color="ACvsW6"), alpha=0.5) +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                              group=SCB_C, color="ACvsW6"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                              group=SCB_C, color="ACvsW6"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS,
                              group=SCB_C, color="ACvsW6"), alpha=0.5) +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                              group=SCB_C, color="ACvsW6"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                              group=SCB_C, color="ACvsW6"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS,
                              group=SCB_C, color="ACvsW1"), alpha=0.5) +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                              group=SCB_C, color="ACvsW1"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                              group=SCB_C, color="ACvsW1"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS,
                              group=SCB_C, color="ACvsD10"), alpha=0.5) +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                              group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                              group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS,
                              group=SCB_C, color="ACvsD10"), alpha=0.5) +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                              group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                              group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS,
                              group=SCB_C, color="ACvsD10"), alpha=0.5) +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                              group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                          geom_line(data=mscACvsrest_AW.ag[SCA_C=="A" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                              group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                                geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS,
                                    group=SCB_C, color="ACvsW6"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="ACvsW6"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.1"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="ACvsW6"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS,
                                      group=SCB_C, color="ACvsW6"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="ACvsW6"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W6_6.3"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="ACvsW6"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS,
                                      group=SCB_C, color="ACvsW1"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="ACvsW1"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Spring_2016_W1_1.1"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="ACvsW1"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS,
                                      group=SCB_C, color="ACvsD10"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="D"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS,
                                      group=SCB_C, color="ACvsD10"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="AF"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS,
                                      group=SCB_C, color="ACvsD10"), alpha=0.5) +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS,
                                      group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                                  geom_line(data=mscACvsrest_AW.ag[SCA_C=="C" & SCB_C=="OO_Fall_2016_D10_41"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS,
                                      group=SCB_C, color="ACvsD10"), alpha=0.5, linetype="dotted") +
                                      geom_line(data=mscACvsrest_AW.ag[SCB_C=="C"], aes(x=as.factor(window), y=meanIBS, group=SCB_C, color="AvsC"), size=1) +
                                          geom_line(data=mscACvsrest_AW.ag[SCB_C=="C"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS, group=SCB_C,
                                            color="AvsC"), linetype="dotted") +
                                          geom_line(data=mscACvsrest_AW.ag[SCB_C=="C"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS, group=SCB_C,
                                            color="AvsC"), linetype="dotted")



    mscACvsrest_A.ag <- mscACvsrest_A[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
      list(window, chr, start, stop, SCcomp, speciesA, speciesB, popA, popB, SCA_C, SCB_C)]

    mscACvsrest_A.agW <- mscACvsrest_A.ag[popA=="W1" | popA=="W6" | popB=="W1" | popB=="W6"]
    mscACvsrest_A.agW$popcomp <- paste(mscACvsrest_A.agW$popA, mscACvsrest_A.agW$popB, sep="_")

  # First thing we want to do is compare AvsC to AvsW and CvsW, or AvsD10 and CvsD10, etc. Let's start with A and C vs W1 and W6
  # Let's make a file with just the W ponds and A and C clones, and remove DMud
    mscACvsrest.agW <- mscACvsrest.ag[popA=="W1" | popA=="W6" | popB=="W1" | popB=="W6"]
    mscACvsrest.agW <- mscACvsrest.agW[popA!="DMud" & popB!="DMud"]

    save(msc, file="m_IBSbyslidingwindow_subsetSC_1000000_100000_nosub_202003025_withsc.Rdata")

  mscwithinsc <- msc[SCA==SCB]
  save(mscwithinsc, file="m_IBSbyslidingwindow_subsetSC_1000000_100000_nosub_202003025_withsc_withinsc.Rdata")
  mscwithinscOO <- mscwithinsc[SCA=="OO"]
  save(mscwithinscOO, file="m_IBSbyslidingwindow_subsetSC_1000000_100000_nosub_202003025_withsc_OOcomparisons.Rdata")




  scsubA <- data.table(cloneA=sc$clone, medrdA=sc$medrd)
  scsubB <- data.table(cloneB=sc$clone, medrdB=sc$medrd)
  setkey(scsubB, cloneB)
  setkey(mscwithinsc, cloneB)
  m <- merge(mscwithinsc, scsubB)
  setkey(scsubA, cloneA)
  setkey(m, cloneA)
  m2 <- merge(m, scsubA)
  m2noOO <- m2[SCA!="OO"]

  m2noOO.ag <- m2noOO[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
    list(SCA, window, chr, start, stop)]

  m2noOO$clonecompare <- paste(m2noOO$cloneA,m2noOO$cloneB,sep="_")
  m2noOO.agbyclone <- m2noOO[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
    list(clonecompare, SCA, medrdA, medrdB)]





  #!/usr/bin/env Rscript

  ### libraries
    library(data.table)
    library(foreach)
    library(ggplot2)

  load("m_IBSbyslidingwindow_subsetSC_1000000_100000_nosub_202003025_withsc.Rdata")
  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexandObtusaC/CloneInfoFilePulexandObtusa_withmedrd_update20200324")

  mscbetweenpop <- msc[SCA!=SCB]

  mscbetweenpop$populationA <- ifelse(mscbetweenpop$populationA=="Dcat", "DCat", mscbetweenpop$populationA)
  mscbetweenpop$populationB <- ifelse(mscbetweenpop$populationB=="Dcat", "DCat", mscbetweenpop$populationB)

  mscbetweenpop$SCcompare <- paste(mscbetweenpop$SCA,mscbetweenpop$SCB, sep="_")
  mscbetweenpop$popcompare <- paste(mscbetweenpop$populationA,mscbetweenpop$populationB, sep="_")
  mscbetweenpop$SCpop <- paste(mscbetweenpop$SCcompare,mscbetweenpop$popcompare,sep="_")

  msc.ag <- mscbetweenpop[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
    list(SCpop, window, chr, start, stop)]

  save(msc.ag, file="m_IBSbyslidingwindow_subsetSC_1000000_100000_nosub_202003025_withsc_betweensc_agg")

  temp <- unlist(strsplit(as.character(msc.ag$SCpop), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  msc.ag$SCcompare <- paste(matdat$V1,matdat$V2,sep="_")
  msc.ag$populationA <- matdat$V3
  msc.ag$populationB <- matdat$V4
  msc.ag$pondcompare <- paste(matdat$V3,matdat$V4,sep="_")

  msc.ag$pondcompare <- ifelse(msc.ag$pondcompare=="D8_D10", "D10_D8", msc.ag$pondcompare)
  msc.ag$pondcompare <- ifelse(msc.ag$pondcompare=="DBunk_D10", "D10_DBunk", msc.ag$pondcompare)
  msc.ag$pondcompare <- ifelse(msc.ag$pondcompare=="DCat_D10", "D10_DCat", msc.ag$pondcompare)
  msc.ag$pondcompare <- ifelse(msc.ag$pondcompare=="DBunk_D8", "D8_DBunk", msc.ag$pondcompare)
  msc.ag$pondcompare <- ifelse(msc.ag$pondcompare=="DCat_D8", "D8_DCat", msc.ag$pondcompare)
  msc.ag$pondcompare <- ifelse(msc.ag$pondcompare=="W1_D8", "D8_W1", msc.ag$pondcompare)
  msc.ag$pondcompare <- ifelse(msc.ag$pondcompare=="W6_D8", "D8_W6", msc.ag$pondcompare)
  msc.ag$pondcompare <- ifelse(msc.ag$pondcompare=="W1_DBunk", "DBunk_W1", msc.ag$pondcompare)
  msc.ag$pondcompare <- ifelse(msc.ag$pondcompare=="W6_DBunk", "DBunk_W6", msc.ag$pondcompare)

  msc.agsub <- msc.ag[populationA!="DOil" & populationA!="Dramp" & populationB!="DOil" & populationB!="Dramp" &
    populationA!="DMud" & populationA!="DLily" & populationB!="DMud" & populationB!="DLily"]

  msc.agsubAC <- msc.agsub[SCcompare=="A_C" | SCcompare=="C_A"]
  msc.agsubnoAC <- msc.agsub[SCcompare!="A_C" | SCcompare!="C_A"]

  msc.agponds <- msc.agsubnoAC[,list(meanIBS=mean(meanIBS), sd_IBS=sd(meanIBS)),
    list(pondcompare, window, chr, start, stop)]

  msc.agAC <- msc.agsubAC[,list(meanIBS=mean(meanIBS), sd_IBS=sd(sd_IBS)),
    list(window, chr, start, stop)]
  msc.agAC$pondcompare <- c("A_C")

  msc.ag2 <- rbind(msc.agAC, msc.agponds)

  msc.ag2$category <- ifelse(msc.ag2$pondcompare=="D8_DBunk" | msc.ag2$pondcompare=="D8_DCat" |
    msc.ag2$pondcompare=="DBunk_DCat", "D8surround", ifelse(msc.ag2$pondcompare=="D10_D8" |
    msc.ag2$pondcompare=="D10_DBunk" | msc.ag2$pondcompare=="D10_DCat", "D10vsD8Surround", ifelse(
    msc.ag2$pondcompare=="D8_W1" | msc.ag2$pondcompare=="D8_W6" |
    msc.ag2$pondcompare=="DBunk_W1" | msc.ag2$pondcompare=="DBunk_W6" |
    msc.ag2$pondcompare=="W1_DCat" | msc.ag2$pondcompare=="W6_DCat", "WvsD8Surround", ifelse(
    msc.ag2$pondcompare=="A_C", "zAvsB", ifelse(msc.ag2$pondcompare=="W1_D10" |
    msc.ag2$pondcompare=="W6_D10" | msc.ag2$pondcompare=="W6_W1", 6, "WithinPond")))))

    ggplot(data=msc.ag2[category!="6"], aes(x=as.factor(window), y=meanIBS, group=as.factor(pondcompare), color=as.factor(category))) + geom_line()
    ggplot(data=msc.ag2[category!="6" & category!="zAvsB"], aes(x=as.factor(window), y=sd_IBS, group=as.factor(pondcompare), color=as.factor(category))) + geom_line()

  msc.ag2$UCI <- msc.ag2$meanIBS+(2*msc.ag2$sd_IBS)

    ggplot() + geom_line(data=msc.ag2[pondcompare=="D8_D8"], aes(x=as.factor(window), y=meanIBS, group=as.factor(pondcompare))) +
      geom_line(data=msc.ag2[pondcompare=="D8_D8"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS, group=as.factor(pondcompare), color="grey")) +
      geom_line(data=msc.ag2[pondcompare=="D8_D8"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS, group=as.factor(pondcompare), color="grey")) +
      geom_line(data=msc.ag2[pondcompare=="A_C"], aes(x=as.factor(window), y=meanIBS, group=as.factor(pondcompare), color="blue"))

      ggplot() + geom_line(data=msc.ag2[pondcompare=="D8_DBunk"], aes(x=as.factor(window), y=meanIBS, group=as.factor(pondcompare))) +
        geom_line(data=msc.ag2[pondcompare=="D8_DBunk"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS, group=as.factor(pondcompare), color="grey")) +
        geom_line(data=msc.ag2[pondcompare=="D8_DBunk"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS, group=as.factor(pondcompare), color="grey")) +
        geom_line(data=msc.ag2[pondcompare=="A_C"], aes(x=as.factor(window), y=meanIBS, group=as.factor(pondcompare), color="blue"))

        ggplot() + geom_line(data=msc.ag2[pondcompare=="D10_D8"], aes(x=as.factor(window), y=meanIBS, group=as.factor(pondcompare))) +
          geom_line(data=msc.ag2[pondcompare=="D10_D8"], aes(x=as.factor(window), y=meanIBS+2*sd_IBS, group=as.factor(pondcompare), color="grey")) +
          geom_line(data=msc.ag2[pondcompare=="D10_D8"], aes(x=as.factor(window), y=meanIBS-2*sd_IBS, group=as.factor(pondcompare), color="grey")) +
          geom_line(data=msc.ag2[pondcompare=="A_C"], aes(x=as.factor(window), y=meanIBS, group=as.factor(pondcompare), color="blue"))


  #msc <- msc[populationA!="DOil" & populationA!="Dramp" & populationB!="DOil" & populationB!="Dramp"]

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
