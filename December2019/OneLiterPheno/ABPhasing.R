### libraries
      library(data.table)
      library(ggplot2)
      library(foreach)
      library(lattice)
      library(tidyr)
      library(tidyverse)
      library(gdsfmt)
      library(SNPRelate)
      library(SeqArray)

### Load parsed file

    phase <- fread("phasingworkingB")
    colnames(phase) <- c("chr", "pos", "ref", "alt", "A", "C")

    phase$Aa <- sapply(strsplit(phase$A,":"), `[`,1)
    phase$Cc <- sapply(strsplit(phase$C,":"), `[`,1)
    phase <- phase[, c("chr", "pos", "ref", "alt", "Aa", "Cc"), with=FALSE]

    phase$A1 <- sapply(strsplit(phase$Aa, split=NULL), `[`,1)
    phase$A2 <- sapply(strsplit(phase$Aa, split=NULL), `[`,3)
    phase$Aphase <- sapply(strsplit(phase$Aa, split=NULL), `[`,2)

    phase$C1 <- sapply(strsplit(phase$Cc, split=NULL), `[`,1)
    phase$C2 <- sapply(strsplit(phase$Cc, split=NULL), `[`,3)
    phase$Cphase <- sapply(strsplit(phase$Cc, split=NULL), `[`,2)

    phase$Atype <- ifelse(phase$A1==phase$A2, "hom", ifelse(phase$A1!=phase$A2 &
      phase$Aphase=="|", "phasehet", "unphasehet"))

    phase$Ctype <- ifelse(phase$C1==phase$C2, "hom", ifelse(phase$C1!=phase$C2 &
      phase$Cphase=="|", "phasehet", "unphasehet"))

    save(phase, file="phase_AandCworking_20200706.Rdata")

### Load in the vcf and superclone files
    genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann.seq.gds")
    sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

    snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
        chr = seqGetData(genofile, "chromosome"),
        pos = seqGetData(genofile, "position"))

    setkey(snps, chr, pos)

    chr.ag <- snps[,list(start = min(pos), stop=max(pos)), list(chr)]
    chr.ag$length <- chr.ag$stop-chr.ag$start
    chrtouse <- chr.ag$chr[chr.ag$length>1000000]

    window.size <- 100000
    step.size <- 100000
    wins <- foreach(chr.i=chrtouse, .combine="rbind")%dopar%{
      #chr.i <- "Scaffold_1863_HRSCAF_2081"
      snp.dt.ag <- snps[J(chr.i)][,list(start=min(pos), stop=max(pos))]

      data.table(chr=chr.i, start=seq(from=snp.dt.ag$start,
                                      to=snp.dt.ag$stop - window.size,
                                      by=step.size),
                            stop =seq(from=snp.dt.ag$start,
                                      to=snp.dt.ag$stop - window.size,
                                      by=step.size) + window.size)
    }

    Cphase <- phase[, c("chr", "pos", "ref", "alt", "C1", "C2", "Cphase", "Ctype"), with=FALSE]
    setkey(Cphase, chr, pos)
    Cphase <- Cphase[Ctype=="phasehet"]

    # First subset used 10 per window, let's try 5 per window as a second go

    Csnpssub <- foreach(win.i =1:(dim(wins)[1]), .errorhandling="remove", .combine="rbind")%dopar%{

      print(paste("Getting data: ", win.i, sep=""))

      snp.dt.tmp <- Cphase[J(wins[win.i]$chr)][pos>=wins[win.i]$start & pos<=wins[win.i]$stop]

      snps.sub <- snp.dt.tmp[sample(nrow(snp.dt.tmp), 5), ]

      snps.sub

      }

### So first let's work on just the C clonal lineage, which we will use for the selfed Cs.

    Cphasesnps <- Csnpssub[, c("chr", "pos")]
    setkey(Cphasesnps, chr, pos)
    setkey(snps, chr, pos)
    Csnps <- merge(Cphasesnps, snps)
    Csnpsids <- Csnps$variant.ids

  # For now let's keep only phased heterozygous sites, keep in mind that in the phased set 1 means an alt allele, while a 2 in a dosage call means 2 reference


  # Pull out the selfed C lab generated clones from the sc file
    selfedC <- sc[SC=="selfedC"]

  # For now let's remove the individual with a medrd of 1
    selfedC <- selfedC[medrd >1]
    setkey(selfedC, clone)
    selfedCids <- selfedC$clone

  # Filter genotype file
    seqSetFilter(genofile, variant.id=Csnpsids)
    seqSetFilter(genofile, sample.id=selfedCids)

  # Pull out genotypes (dosage)
    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snps, variant.ids)

    mhet <- merge(snps, het)

    mhetlong <- melt(mhet, measure.vars=selfedCids, variable.name="clone", value.name="dosage")

    setkey(mhetlong, chr, pos)
    setkey(Cphase, chr, pos)
    Cmhetlong <- merge(mhetlong, Cphase)

    Cmhetlong$rqtldosage <- ifelse(Cmhetlong$dosage==1, "H", ifelse(Cmhetlong$dosage==0 &
      Cmhetlong$C1==1, "A", ifelse(Cmhetlong$dosage==2 & Cmhetlong$C1==0, "A", ifelse(Cmhetlong$dosage==0 &
      Cmhetlong$C2==1, "B", ifelse(Cmhetlong$dosage==2 & Cmhetlong$C2==0, "B", "-")))))

    Cmhetlong$snpid <- paste(Cmhetlong$chr, Cmhetlong$pos, sep="_")
    Cmhetlongsubids <- Cmhetlong[, c("snpid", "chr", "pos")]
    Cmhetlongsubids <- unique(Cmhetlongsubids)

    Cmhetlongsub <- Cmhetlong[, c("snpid", "chr", "pos", "clone", "rqtldosage")]

    Cphasewide <- dcast(Cmhetlongsub, snpid ~ clone, value.var="rqtldosage")
    setkey(Cphasewide, snpid)
    setkey(Cmhetlongsubids, snpid)
    mCphasewise <- merge(Cmhetlongsubids, Cphasewide)

    mCphasewisem <- as.matrix(mCphasewise)
    mCphasewisemt <- t(mCphasewisem)
    mCphasewisemtdt <- as.data.table(mCphasewisemt)

    write.csv(mCphasewise, file="selfedCforrqtlsub_5.csv", row.names=FALSE, quote=FALSE)

### Ok, now let's try the AxC F1 hybrids, for now let's remove all sites that are unphased hets in A or C
  # Also remove all sites that are fixed homozygous between A and C
    ACphase <- phase[Atype!="unphasehet" & Ctype!="unphasehet"]
    ACphase$fixhom <- ifelse(ACphase$Atype=="hom" & ACphase$Ctype=="hom", 1, 0)
    ACphase <- ACphase[fixhom!=1]
    ACphase <- ACphase[, c("chr", "pos", "ref", "alt", "A1", "A2", "Aphase", "Atype", "C1", "C2",
      "Cphase", "Ctype"), with=FALSE]

    setkey(ACphase, chr, pos)

    ACsnpssub <- foreach(win.i =1:(dim(wins)[1]), .errorhandling="remove", .combine="rbind")%dopar%{

      print(paste("Getting data: ", win.i, sep=""))

      snp.dt.tmp <- ACphase[J(wins[win.i]$chr)][pos>=wins[win.i]$start & pos<=wins[win.i]$stop]

      snps.sub <- snp.dt.tmp[sample(nrow(snp.dt.tmp), 10), ]

      snps.sub

      }


    ACphasesnps <- ACsnpssub[, c("chr", "pos")]
    setkey(ACphasesnps, chr, pos)
    setkey(snps, chr, pos)
    ACsnps <- merge(ACphasesnps, snps)
    ACsnpsids <- ACsnps$variant.ids

  # Pull out the AxB F1 hybrids that are phenotyped from the sc file
    ACF1 <- sc[OneLiterPheno==1 & AxCF1Hybrid==1]

  # Ok, how to deal with clonal lineages that where phenotyped/genotyped more than once
  # For now, let's take the individual with the highest read depth for S and AL

    ACF1Ss <- ACF1[SC=="S"]
    ACF1ALs <- ACF1[SC=="AL"]
    ACF1other <- ACF1[SC!="AL" & SC!="S"]
    ACF1Sssub <- ACF1Ss[medrd>11]
    ACF1ALssub <- ACF1ALs[medrd>5]
    ACF1subfull <- rbind(ACF1Sssub, ACF1ALssub, ACF1other)
    ACF1subfullids <- ACF1subfull$clone


  # Filter genotype file
    seqResetFilter(genofile)

    seqSetFilter(genofile, variant.id=ACsnpsids)
    seqSetFilter(genofile, sample.id=ACF1subfullids)

  # Pull out genotypes (dosage)
    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(snps, variant.ids)

    mhet <- merge(snps, het)

    mhetlong <- melt(mhet, measure.vars=ACF1subfullids, variable.name="clone", value.name="dosage")

    setkey(mhetlong, chr, pos)
    setkey(ACphase, chr, pos)
    ACmhetlong <- merge(mhetlong, ACphase)

    ACmhetlong$rqtldosage <- ifelse(ACmhetlong$dosage==1 & ACmhetlong$A1==0 & ACmhetlong$A2==1 &
      ACmhetlong$C1==0 & ACmhetlong$C2==0, 6, ifelse(ACmhetlong$dosage==1 & ACmhetlong$A1==1 &
      ACmhetlong$A2==0 & ACmhetlong$C1==0 & ACmhetlong$C2==0, 5, ifelse(ACmhetlong$dosage==1 &
      ACmhetlong$A1==0 & ACmhetlong$A2==0 & ACmhetlong$C1==1 & ACmhetlong$C2==0, 7, ifelse(
      ACmhetlong$dosage==1 & ACmhetlong$A1==0 & ACmhetlong$A2==0 & ACmhetlong$C1==0 &
      ACmhetlong$C2==1, 8, 0))))

    ACmhetlong$rqtldosage <- ifelse(ACmhetlong$dosage==0 & ACmhetlong$A1==0 & ACmhetlong$A2==1 &
      ACmhetlong$C1==0 & ACmhetlong$C2==0, 5, ifelse(ACmhetlong$dosage==0 & ACmhetlong$A1==1 &
      ACmhetlong$A2==0 & ACmhetlong$C1==0 & ACmhetlong$C2==0, 6, ifelse(ACmhetlong$dosage==0 &
      ACmhetlong$A1==0 & ACmhetlong$A2==0 & ACmhetlong$C1==1 & ACmhetlong$C2==0, 8, ifelse(
      ACmhetlong$dosage==0 & ACmhetlong$A1==0 & ACmhetlong$A2==0 & ACmhetlong$C1==0 &
      ACmhetlong$C2==1, 7, ACmhetlong$rqtldosage
      ))))

    ACmhetlong$rqtldosage <- ifelse(ACmhetlong$dosage==2 & ACmhetlong$A1==0 & ACmhetlong$A2==1 &
      ACmhetlong$C1==0 & ACmhetlong$C2==0, 5, ifelse(ACmhetlong$dosage==2 & ACmhetlong$A1==1 &
      ACmhetlong$A2==0 & ACmhetlong$C1==0 & ACmhetlong$C2==0, 6, ifelse(ACmhetlong$dosage==2 &
      ACmhetlong$A1==0 & ACmhetlong$A2==0 & ACmhetlong$C1==1 & ACmhetlong$C2==0, 8, ifelse(
      ACmhetlong$dosage==2 & ACmhetlong$A1==0 & ACmhetlong$A2==0 & ACmhetlong$C1==0 &
      ACmhetlong$C2==1, 7, ACmhetlong$rqtldosage
      ))))

    ACmhetlong$snpid <- paste(ACmhetlong$chr, ACmhetlong$pos, sep="_")
    ACmhetlongsubids <- ACmhetlong[, c("snpid", "chr", "pos")]
    ACmhetlongsubids <- unique(ACmhetlongsubids)

    ACmhetlongsub <- ACmhetlong[, c("snpid", "chr", "pos", "clone", "rqtldosage")]

    ACphasewide <- dcast(ACmhetlongsub, snpid ~ clone, value.var="rqtldosage")
    setkey(ACphasewide, snpid)
    setkey(ACmhetlongsubids, snpid)
    mACphasewise <- merge(ACmhetlongsubids, ACphasewide)

    mACphasewisem <- as.matrix(mACphasewise)
    mACphasewisemt <- t(mACphasewisem)
    mACphasewisemtdt <- as.data.table(mACphasewisemt)

    write.csv(mACphasewise, file="AxCF1forrqtlsub.csv", row.names=FALSE, quote=FALSE)
    write.table(mACphasewisemtdt, file="AxCF1forrqtlt.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)


    A <- fread("AxCF1forrqtlTA.csv")
    B <- fread("AxCF1forrqtlTB.csv")
    C <- fread("AxCF1forrqtlTC.csv")
    AB <- cbind(A, B, C)
    AB[is.na(AB)] <- "-"

    write.csv(AB, file="AxCF1forrqtltake2.csv", row.names=FALSE, quote=FALSE)




    ACF1Ss <- ACF1[SC=="S"]
    ACF1Ssids <- ACF1Ss$clone

  # Filter genotype file
    seqSetFilter(genofile, variant.id=ACsnpsids)
    seqSetFilter(genofile, sample.id=ACF1Ssids)

  # Pull out genotypes (dosage)
    het <- t(seqGetData(genofile, "$dosage"))
    het <- as.data.table(het)

    colnames(het) <- c(seqGetData(genofile, "sample.id"))
    het$variant.ids <- seqGetData(genofile, "variant.id")

    setkey(het, variant.ids)
    setkey(ACsnps, variant.ids)

    mhet <- merge(ACsnps, het)

    mhetlong <- melt(mhet, measure.vars=ACF1Ssids, variable.name="clone", value.name="dosage")

    clonecounts <- mhetlong[, .N, by=list(variant.ids, dosage)]
    clonecounts <- clonecounts[dosage!="NA"]

    clonecountswide <- dcast(clonecounts, variant.ids ~ dosage, value.var="N")
    colnames(clonecountswide) <- c("clone", "dos0", "dos1", "dos2")
    clonecountswide[is.na(dos0),dos0:=0]
    clonecountswide[is.na(dos1),dos1:=0]
    clonecountswide[is.na(dos2),dos2:=0]
    clonecountswide$total <- clonecountswide$dos0 + clonecountswide$dos1 + clonecountswide$dos2
    clonecountswide$consensusS <- ifelse(clonecountswide$dos0 > clonecountswide$dos1 & clonecountswide$dos0 >
      clonecountswide$dos2, 0, ifelse(clonecountswide$dos1>clonecountswide$dos0 & clonecountswide$dos1 >
      clonecountswide$dos2, 1, ifelse(clonecountswide$dos2 > clonecountswide$dos0 & clonecountswide$dos2 >
      clonecountswide$dos1, 2, "other")))

      ### libraries
            library(data.table)
            library(ggplot2)
            library(foreach)
            library(lattice)
            library(tidyr)
            library(tidyverse)
            library(gdsfmt)
            library(SNPRelate)
            library(SeqArray)

    phaserc <- fread("AxC_F1.csv")

    phaserc <- phaserc[, c("clone", "chr.x", "pos", "founder")]
    colnames(phaserc) <- c("clone", "chr", "pos", "founder")

    phasercsnps <- phaserc[, c("chr", "pos")]
    phasercsnps <- unique(phasercsnps)

    setkey(phasercsnps, chr, pos)

  ### Load in the vcf and superclone files
    genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann.seq.gds")

    snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
        chr = seqGetData(genofile, "chromosome"),
        pos = seqGetData(genofile, "position"))

    setkey(snps, chr, pos)

    chr.ag <- snps[,list(start = min(pos), stop=max(pos)), list(chr)]
    chr.ag$length <- chr.ag$stop-chr.ag$start
    chrtouse <- chr.ag$chr[chr.ag$length>1000000]

    window.size <- 100000
    step.size <- 100000
    wins <- foreach(chr.i=chrtouse, .combine="rbind")%dopar%{
      #chr.i <- "Scaffold_1863_HRSCAF_2081"
      snp.dt.ag <- snps[J(chr.i)][,list(start=min(pos), stop=max(pos))]

      data.table(chr=chr.i, start=seq(from=snp.dt.ag$start,
                                      to=snp.dt.ag$stop - window.size,
                                      by=step.size),
                            stop =seq(from=snp.dt.ag$start,
                                      to=snp.dt.ag$stop - window.size,
                                      by=step.size) + window.size)
    }


  # Try with 3 SNPs per window

    ACrcsub <- foreach(win.i =1:(dim(wins)[1]), .errorhandling="remove", .combine="rbind")%dopar%{

      print(paste("Getting data: ", win.i, sep=""))

      snp.dt.tmp <- phasercsnps[J(wins[win.i]$chr)][pos>=wins[win.i]$start & pos<=wins[win.i]$stop]

      snps.sub <- snp.dt.tmp[sample(nrow(snp.dt.tmp), 3), ]

      snps.sub

      }

    setkey(ACrcsub, chr, pos)
    setkey(phaserc, chr, pos)
    subphaserc <- merge(ACrcsub, phaserc)

    subphaserc$snpid <- paste(subphaserc$chr,subphaserc$pos,sep="_")
    subphaserc$rqtlgeno <- ifelse(subphaserc$founder=="A", 1, ifelse(subphaserc$founder=="B", 2, ifelse(
      subphaserc$founder=="C", 3, 4
    )))

    subphaserc <- subphaserc[, c("chr", "pos", "snpid", "clone", "rqtlgeno"), with=FALSE]
    subphasercids <- subphaserc[, c("snpid", "chr", "pos")]
    subphasercids <- unique(subphasercids)

    phasercwide <- dcast(subphaserc, snpid ~ clone, value.var="rqtlgeno")
    setkey(phasercwide, snpid)
    setkey(subphasercids, snpid)
    mCphasewise <- merge(subphasercids, phasercwide)
    setkey(mCphasewise, chr, pos)

    write.csv(mCphasewise, file="AxCF1recontruct_sub3.csv", row.names=FALSE, quote=FALSE)
