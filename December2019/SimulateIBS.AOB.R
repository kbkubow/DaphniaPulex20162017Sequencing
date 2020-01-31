
#ijob -c1 -p standard -A berglandlab
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

#!/usr/bin/env Rscript

### libraries
  library(gdsfmt)
  library(SNPRelate)
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(lattice)
  library(tidyr)
  library(SeqArray)

### set working directory
  setwd("/project/berglandlab/Karen/MappingDec2019/")

### open GDS file
  genofile <- seqOpen("MapDec19PulexOnlyB_filtsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### Load SuperClone file
  sc <- fread("CloneInfoFilePulexOnly_20200128.txt")

### Load LD-pruned SNP-set
  load("finalsetsnpset01pulex_20200121.Rdata")
  seqSetFilter(genofile, variant.id=finalsetsnpset01)

  snps <- data.table(variant.id = seqGetData(genofile, "variant.id"),
    chr = seqGetData(genofile, "chromosome"),
    pos = seqGetData(genofile, "position"))

### Simulation function
  simGeno <- function(SC.i="B", nSNPs=length(finalsetsnpset01)) {
    #SC.i <- "B"; nSNPs <- 100

    ### first, choose random individual from super clone set
      focal.individual <- sample(sc[SC==SC.i]$clone, 1)

    ### get focal individual's genotype
      seqResetFilter(genofile)
      seqSetFilter(genofile, sample.id=focal.individual, variant.id=finalsetsnpset01[1:nSNPs])

      focal.genotype <- cbind(snps[1:nSNPs], data.table(af=abs(1-seqAlleleFreq(genofile))))

    ### get empirical read depths for all other members of that superclone
      seqResetFilter(genofile)
      seqSetFilter(genofile, sample.id=sc[SC==SC.i][clone!=focal.individual]$clone, variant.id=finalsetsnpset01[1:nSNPs])

      obs.rd <- cbind(snps[1:nSNPs], as.data.frame(t(seqGetData(genofile, "annotation/format/DP")$data)))
      setnames(obs.rd, names(obs.rd)[grepl("V[0-9]*", names(obs.rd))], seqGetData(genofile, "sample.id"))
      obs.rd <- melt(obs.rd, id.vars=c("variant.id", "chr", "pos"), variable.name="clone", value.name="rd")

    ### get empirical genotype
      obs.geno <- cbind(snps[1:nSNPs], as.data.frame(t(seqGetData(genofile, "$dosage"))))
      setnames(obs.geno, names(obs.geno)[grepl("V[0-9]*", names(obs.geno))], seqGetData(genofile, "sample.id"))
      obs.geno <- melt(obs.geno, id.vars=c("variant.id", "chr", "pos"), variable.name="clone", value.name="obsGeno")

    ### determine "true" simulated genotype of focal individual
      setkey(focal.genotype, variant.id, chr, pos)
      setkey(obs.rd, variant.id, chr, pos)
      setkey(obs.geno, variant.id, chr, pos)

      obs.rd <- merge(obs.rd, focal.genotype)

      setkey(obs.geno, variant.id, chr, pos, clone)
      setkey(obs.rd, variant.id, chr, pos, clone)
      obs.rd <- merge(obs.rd, obs.geno)

    ### remove individuals with missing genotype calls the first time arounbd
      obs.rd <- obs.rd[!is.na(af) & !is.na(obsGeno)]

    ### simulate the observed genotype based on empirical read depth & "true" genotype
      obs.rd[,sim.nAlt:=rbinom(dim(obs.rd)[1], obs.rd$rd, obs.rd$af)]
      obs.rd[,sim.nRef:=rd - sim.nAlt]

    ### include GATK look up table at some point
      obs.rd[sim.nAlt==0 & sim.nRef>0, simGeno:=0]
      obs.rd[sim.nAlt>0 & sim.nRef>0, simGeno:=1]
      obs.rd[sim.nAlt>0 & sim.nRef==0, simGeno:=2]

    ### return
      obs.rd
    }

### IBS function
  ibs.fun <- function(dt) {
    #dt <- sim

    uniq.clones <- unique(dt$clone)
    setkey(dt, clone)

    o <- foreach(i=1:(length(uniq.clones) - 1), .combine="rbind")%do%{
      foreach(j=(i+1):length(uniq.clones), .combine="rbind")%do%{
        print(paste(i, j, sep=" / "))

        tmp <- merge(dt[J(uniq.clones[i])], dt[J(uniq.clones[j])], by="variant.id")


        data.table(cloneA=uniq.clones[i],
                   cloneB=uniq.clones[j],
                   sim.IBS=1-(mean(abs(tmp$simGeno.x - tmp$simGeno.y), na.rm=TRUE)/2),
                   obs.IBS=1-(mean(abs(tmp$obsGeno.x - tmp$obsGeno.y), na.rm=TRUE)/2))

      }
    }
  }

### run
  foreach(i=c("A", "B", "H")) {
    sim <- simGeno(SC.i="H")
    ibs.out <- ibs.fun(sim)

  }

  summary(ibs.out)

### save
  save(ibs.out, sim, file="/nv/vol186/bergland-lab/alan/ibs.out.Rdata")



#### on workstatiopn
  library(data.table)
  library(ggplot2)
  library(ggbeeswarm)

  load("/mnt/sammas_storage/bergland-lab/alan/ibs.out.Rdata")

  ibs.out.l <- melt(ibs.out, id.vars=c("cloneA", "cloneB"))

  ggplot(data=ibs.out.l, aes(x=variable, y=value, color=variable)) + geom_beeswarm()


  sfs <- sim[,list(obs.f=mean(abs(2-obsGeno)/2), sim.f=mean(simGeno)/2), list(variant.id)]
  sfsl <- melt(sfs, id.vars="variant.id")
  sfsl[,fold:=value]
  sfsl[value>.5, fold:=1-value]

  ggplot(data=sfsl[fold>.4], aes(fold, group=variable)) + geom_histogram() + facet_wrap(~variable)

  wilcox.test(fold~variable, sfsl[fold>.4])


  sfsw <- dcast(sfsl, variant.id~variable, value.var="fold")


# Let's pull out genotypes for all of these
    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snps, variant.ids)

    mhet <- merge(snps, het)

    mhetlong <- melt(mhet, measure.vars=scBids, variable.name="clone", value.name="dosage")

# Now let's pull out the read depths
    dp <- t((seqGetData(genofile, "annotation/format/DP"))$data)
    dp <- as.data.table(dp)

    colnames(dp) <- c(seqGetData(genofile, "sample.id"))
    dp$variant.ids <- seqGetData(genofile, "variant.id")

    dplong <- melt(dp, measure.vars=scBids, variable.name="clone", value.name="dp")

# Now let's choose a starting reference clone and load the gatk probability table
      mhetlongref <- mhetlong[clone=="April_2017_D8_125"]
      mhetlongrefsub <- mhetlongref[, c("variant.ids", "dosage"), with=FALSE]
      setkey(mhetlongrefsub, variant.ids)

      gatkprob <- fread("gtProbabilities.txt")
      gtprobfilt <- gatkprob[2:3212]
      gtprobfiltsub <- gtprobfilt[,c("ref.rd", "alt.rd", "GTprob2", "GTprob1", "GTprob0")]


      fullsummarymodIBS <- foreach(x=1:(length(scBids)), .combine="rbind")%do%{
        #x<-1
        testclone <- scBids[x]

        mhetlongrefsub$testrd <- dp[,testclone, with=F]
        mhetlongrefsub$testdosage <- mhet[,testclone, with=F]
        mhetlongrefsub$rbinommod <- rbinom(150408, mhetlongrefsub$testrd, 0.5)
        mhetlongrefsub$ref.rd <- as.integer(ifelse(mhetlongrefsub$dosage==2, mhetlongrefsub$testrd,
          ifelse(mhetlongrefsub$dosage==1, mhetlongrefsub$rbinommod,
          ifelse(mhetlongrefsub$dosage==0, 0, "NA"))))
        mhetlongrefsub$alt.rd <- mhetlongrefsub$testrd-mhetlongrefsub$ref.rd
        mhetlongrefsubB <- mhetlongrefsub[mhetlongrefsub$testdosage==0 | mhetlongrefsub$testdosage==2 |
          mhetlongrefsub$testdosage==2]

        setkey(mhetlongrefsubB, ref.rd, alt.rd)
        setkey(gtprobfiltsub, ref.rd, alt.rd)
        mtmp <- merge(mhetlongrefsubB, gtprobfiltsub)
        mtmp[,Ocall:= sample(c("2","1","0"), prob=c(GTprob2, GTprob1, GTprob0), size=1), by=variant.ids]
        mtmpsub <- mtmp[,c("variant.ids", "Ocall")]
        colnames(mtmpsub) <- c("variant.ids", testclone)
        mhetlongrefsubsub <- mhetlongrefsub[,"variant.ids"]
        setkey(mtmpsub, variant.ids)
        setkey(mhetlongrefsubsub, variant.ids)
        mmtmp <- merge(mhetlongrefsubsub, mtmpsub, all.x=TRUE)
        setkey(mhetlongrefsub, variant.ids)
        mhetlongrefsub[[testclone]] <- mmtmp[[testclone]]

        }

        mhetlongsim <- melt(mhetlongrefsub, measure.vars=scBids, variable.name="clone", value.name="dosage")
        mhetlongsimsub <- mhetlongsim[,c("variant.ids", "clone", "dosage.1")]

  ### Now that we have table of simulated genotypes, need to calculate IBS, should also calculate it for original genotypes
  ### First let's try original genotypes
    # First need to get a list of all pairwise combinations of Bs.

      pairwiseids <- foreach(x=1:(length(scBids)), .combine="rbind")%do%{

        testclone <- scBids[x]
        data.table(cloneA=c(testclone), cloneB=scBids)

        }

    #Remove identical comparisons
      pairwiseidsfilt <- pairwiseids[cloneA!=cloneB]
    #Keep single way comparisons
      pairwiseidsfilt$cloneAnum <- as.numeric(as.factor(pairwiseidsfilt$cloneA))
      pairwiseidsfilt$cloneBnum <- as.numeric(as.factor(pairwiseidsfilt$cloneB))
      pairwiseidsfilt$CloneNumComb <- ifelse(pairwiseidsfilt$cloneAnum > pairwiseidsfilt$cloneBnum,
        paste(pairwiseidsfilt$cloneAnum,pairwiseidsfilt$cloneBnum,sep="_"),
        paste(pairwiseidsfilt$cloneBnum,pairwiseidsfilt$cloneAnum,sep="_"))
      setkey(pairwiseidsfilt, CloneNumComb)
      pairwiseidsfiltunique <- unique(pairwiseidsfilt, by="CloneNumComb")
    #Now get back to the original three columns
      pairwiseidsfiltuniqueB <- data.table(cloneA=pairwiseidsfiltunique$cloneA,
        cloneB=pairwiseidsfiltunique$cloneB)

  ### Ok, now let's calculate IBS for all original B genotypes - have to figure out how IBS is calculated

      IBSoriginal <- foreach(x=1:dim(pairwiseidsfiltuniqueB)[1], .combine="rbind")%do%{
        cloneA=pairwiseidsfiltuniqueB$cloneA[[x]]
        cloneB=pairwiseidsfiltuniqueB$cloneB[[x]]

        tmp <- mhet[, c(cloneA, cloneB), with=F]
        tmp$dist <- abs(tmp[[cloneA]]-tmp[[cloneB]])
        data.table(cloneA=cloneA, cloneB=cloneB, IBS=1-(mean(tmp$dist, na.rm=TRUE)/2))

      }

  ### And now let's calculate IBS for all simulated B genotypes - have to figure out how IBS is calculated

      IBSsimulated <- foreach(x=1:dim(pairwiseidsfiltuniqueB)[1], .combine="rbind")%do%{
        cloneA=pairwiseidsfiltuniqueB$cloneA[[x]]
        cloneB=pairwiseidsfiltuniqueB$cloneB[[x]]

        tmp <- mhetlongrefsub[, c(cloneA, cloneB), with=F]
        tmp$dist <- abs(as.integer(tmp[[cloneA]])-as.integer(tmp[[cloneB]]))
        data.table(cloneA=cloneA, cloneB=cloneB, SimIBS=1-(mean(tmp$dist, na.rm=TRUE)/2))

      }

  ### Now let's put the original and simulated IBS distributions in one table and save
      setkey(IBSoriginal, cloneA, cloneB)
      setkey(IBSsimulated, cloneA, cloneB)
      IBSorigsim <- merge(IBSoriginal, IBSsimulated)
      save(IBSorigsim, file="IBSorigsim_Bclones_20200129.Rdata")

  ### Also added in the original snpgds IBS distribution, saving below
      save(totalIBSwithprogram, file="totalIBSwithprogram_20200129.Rdata")


  ### Now let's look at the allele site frequency spectra for original and simulated genotypes
    # Start with original
      dosagecountsoriginal <- mhetlong[, .N, by=list(variant.ids, dosage)]
      dosagecountsoriginalnoNA <- dosagecountsoriginal[dosage==0 | dosage==1 | dosage==2]

    #Transform to wide format
        dosagecountsoriginalwide <- dcast(dosagecountsoriginalnoNA, variant.ids ~ dosage, value.var="N")
        colnames(dosagecountsoriginalwide) <- c("variant.ids", "dos0", "dos1", "dos2")
        dosagecountsoriginalwide[is.na(dos0),dos0:=0]
        dosagecountsoriginalwide[is.na(dos1),dos1:=0]
        dosagecountsoriginalwide[is.na(dos2),dos2:=0]
        dosagecountsoriginalwide$total <- dosagecountsoriginalwide$dos0+dosagecountsoriginalwide$dos1+
          dosagecountsoriginalwide$dos2
        dosagecountsoriginalwide$propalt <- (dosagecountsoriginalwide$dos0*2 +
          dosagecountsoriginalwide$dos1)/(dosagecountsoriginalwide$total*2)
        dosagecountsoriginalwide$foldedminorallele <- ifelse(dosagecountsoriginalwide$propalt > 0.5,
          1-dosagecountsoriginalwide$propalt, dosagecountsoriginalwide$propalt)

        save(dosagecountsoriginalwide, file="dosagecountsoriginalwide_SCB_20200129.Rdata")

    # Now do the same with simulated
        colnames(mhetlongsimsub) <- c("variant.ids", "clone", "dosage")
        dosagecountssim <- mhetlongsimsub[, .N, by=list(variant.ids, dosage)]
        dosagecountssimnoNA <- dosagecountssim[dosage==0 | dosage==1 | dosage==2]

    #Transform to wide format
        dosagecountssimwide <- dcast(dosagecountssimnoNA, variant.ids ~ dosage, value.var="N")
        colnames(dosagecountssimwide) <- c("variant.ids", "dos0", "dos1", "dos2")
        dosagecountssimwide[is.na(dos0),dos0:=0]
        dosagecountssimwide[is.na(dos1),dos1:=0]
        dosagecountssimwide[is.na(dos2),dos2:=0]
        dosagecountssimwide$total <- dosagecountssimwide$dos0+dosagecountssimwide$dos1+
            dosagecountssimwide$dos2
        dosagecountssimwide$propalt <- (dosagecountssimwide$dos0*2 +
            dosagecountssimwide$dos1)/(dosagecountssimwide$total*2)
        dosagecountssimwide$foldedminorallele <- ifelse(dosagecountssimwide$propalt > 0.5,
            1-dosagecountssimwide$propalt, dosagecountssimwide$propalt)

        save(dosagecountssimwide, file="dosagecountssimwide_SCB_20200129.Rdata")

    # Combine two files together
        dosagecountsoriginalwidesub <- data.table(variant.ids=dosagecountsoriginalwide$variant.ids,
          foldedminorallele=dosagecountsoriginalwide$foldedminorallele, type=c("original"))
        dosagecountssimwidesub <- data.table(variant.ids=dosagecountssimwide$variant.ids,
          foldedminorallele=dosagecountssimwide$foldedminorallele, type=c("simulated"))
        allelespectra <- rbind(dosagecountsoriginalwidesub, dosagecountssimwidesub)

        ggplot(data=allelespectra, aes(x=foldedminorallele, fill=type)) + geom_histogram() + facet_wrap(~type)
        ggplot(data=allelespectra[foldedminorallele!=0 & foldedminorallele!=0.5],
          aes(x=foldedminorallele, fill=type)) + geom_histogram() + facet_wrap(~type)
