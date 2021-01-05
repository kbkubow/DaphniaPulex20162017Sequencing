module load  gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3

R

#!/usr/bin/env Rscript

### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(tidyverse)
  library(gdsfmt)
  library(SNPRelate)
  library(SeqArray)

  genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann.seq.gds")

  filtsnp <- fread("snpsvarpulexpresentinhalf_table_20200623")
  filtsnpids <- filtsnp$variant.ids
  seqSetFilter(genofile, variant.id=filtsnpids)

  snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
      chr = seqGetData(genofile, "chromosome"),
      pos = seqGetData(genofile, "position"),
      dp = seqGetData(genofile, "annotation/info/DP"))

### Narrow down to SNPs that are heterozygous in A or C.

  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  scAC <- sc[SC=="A" | SC=="C" & Nonindependent==0]
  scACids <- scAC$clone
  seqSetFilter(genofile, sample.id=scACids)

# Pull out genotypes
  het <- t(seqGetData(genofile, "$dosage"))
  het <- as.data.table(het)

  colnames(het) <- c(seqGetData(genofile, "sample.id"))
  het$variant.ids <- seqGetData(genofile, "variant.id")

  setkey(het, variant.ids)
  setkey(snps, variant.ids)

  mhet <- merge(snps, het)

  mhetlong <- melt(mhet, measure.vars=scACids, variable.name="clone", value.name="dosage")

  setkey(scAC, clone)
  setkey(mhetlong, clone)
  mmhetlong <- merge(mhetlong, scAC)

  SCcounts <- mmhetlong[, .N, by=list(SC, variant.ids, dosage)]

#Remove NAs
  SCcounts <- SCcounts[dosage!="NA"]
  SCcountsA <- SCcounts[SC=="A"]
  SCcountsC <- SCcounts[SC=="C"]

  SCcountsAwide <- dcast(SCcountsA, variant.ids ~ dosage, value.var="N")
  colnames(SCcountsAwide) <- c("variant.ids", "Ados0", "Ados1", "Ados2")
  SCcountsAwide[is.na(Ados0),Ados0:=0]
  SCcountsAwide[is.na(Ados1),Ados1:=0]
  SCcountsAwide[is.na(Ados2),Ados2:=0]

  SCcountsCwide <- dcast(SCcountsC, variant.ids ~ dosage, value.var="N")
  colnames(SCcountsCwide) <- c("variant.ids", "Cdos0", "Cdos1", "Cdos2")
  SCcountsCwide[is.na(Cdos0),Cdos0:=0]
  SCcountsCwide[is.na(Cdos1),Cdos1:=0]
  SCcountsCwide[is.na(Cdos2),Cdos2:=0]

  setkey(SCcountsAwide, variant.ids)
  setkey(SCcountsCwide, variant.ids)
  SCcountsAC <- merge(SCcountsAwide, SCcountsCwide)

  SCcountsAC$Atot <- SCcountsAC$Ados0+SCcountsAC$Ados1+SCcountsAC$Ados2
  SCcountsAC$Ctot <- SCcountsAC$Cdos0+SCcountsAC$Cdos1+SCcountsAC$Cdos2

  SCcountsACfilt <- SCcountsAC[Atot>59 & Ctot>19]

  SCcountsACfiltAhetCnot <- SCcountsACfilt[Ados1>=Atot-3 & Cdos1<=Ctot-4]
  SCcountsACfiltAnotChet <- SCcountsACfilt[Cdos1>=Ctot-3 & Ados1<=Atot-4]
  SCcountsACfiltAheChet <- SCcountsACfilt[Cdos1>=Ctot-3 & Ados1>=Atot-3]

  ACpoly <- rbind(SCcountsACfiltAhetCnot, SCcountsACfiltAnotChet, SCcountsACfiltAheChet)
  ACpolyids <- ACpoly$variant.ids
  save(ACpolyids, file="ACpolyids_20201117.Rdata")

  seqSetFilter(genofile, variant.id=ACpolyids)

  snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
      chr = seqGetData(genofile, "chromosome"),
      pos = seqGetData(genofile, "position"),
      dp = seqGetData(genofile, "annotation/info/DP"))


### get annotation
  tmp <- seqGetData(genofile, "annotation/info/ANN")
  len1 <- tmp$length
  len2 <- tmp$data
  snp.dt1 <- data.table(len=rep(len1, times=len1),
                ann=len2,
                id=rep(snps$variant.id, times=len1))

  # Extracting data between the 2nd and third | symbol
      snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
      snp.dt1[,impact:=tstrsplit(snp.dt1$ann,"\\|")[[3]]]
      snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]
  # Collapsing additional annotations to original SNP vector length
      snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","),
        impact=paste(impact, collapse=","), gene=paste(gene, collapse=",")),
                        list(variant.ids=id)]
      snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
      snp.dt1.an[,impact:=tstrsplit(snp.dt1.an$impact,"\\,")[[1]]]
      snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]

### Add annotations to snp table
  setkey(snps, variant.ids)
  setkey(snp.dt1.an, variant.ids)
  snpsann <- merge(snps, snp.dt1.an)

  # Load in peaks file
    #peaks <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/peaks.csv")
    load("gprime_peaks.replicates.250K.05.Rdata")

    peaks$narrowstart <- peaks$posMaxGprime-25000
    peaks$narrowstop <- peaks$posMaxGprime+25000


  # Figure out which snps fall in peaks

  snpsinpeaks <- foreach(peak.i=1:dim(peaks)[1], .combine="rbind")%do%{
    #peak.i <- 1
    c <- peaks$CHROM[[peak.i]]
    strt <- peaks$start[[peak.i]]
    stp <- peaks$end[[peak.i]]

    tmp <- snpsann[chr==c & pos >= strt & pos <= stp]
    tmp$peaknum <- c(peak.i)
    tmp
  }

  snpsinpeaksnarrow <- foreach(peak.i=1:dim(peaks)[1], .combine="rbind")%do%{
    #peak.i <- 1
    c <- peaks$CHROM[[peak.i]]
    strt <- peaks$narrowstart[[peak.i]]
    stp <- peaks$narrowstop[[peak.i]]

    tmp <- snpsann[chr==c & pos >= strt & pos <= stp]
    tmp$peaknum <- c(peak.i)
    tmp
  }

### Load in allele freq for male and female pools
  load("alleleFreqs.replicates.250K.05.Rdata")

  alleleFreqs <- alleleFreqs %>%
    rename(
      chr=CHROM,
      pos=POS
    )

  setkey(snpsinpeaksnarrow, chr, pos)
  setkey(alleleFreqs, chr, pos)
  m <- merge(snpsinpeaksnarrow, alleleFreqs)

  m$sex <- ifelse(m$deltaSNP>0, "male", "female")

  m$ALT_FREQ_LOW <- m$AD_ALT.LOW/m$DP.LOW
  m$ALT_FREQ_HIGH <- m$AD_ALT.HIGH/m$DP.HIGH

  mag <- m[, list(difffreqLOW=abs(diff(ALT_FREQ_LOW)), difffreqHIGH=abs(diff(ALT_FREQ_HIGH))), list(chr, pos, peaknum)]
  fractdiv.ag <- fractdiv[,list(minpropdiv=min(propdiv)), list(chr, peaknum)]

  magmean <- m[, list(meanAltFreqLow=mean(ALT_FREQ_LOW), meanAltFreqHigh=mean(ALT_FREQ_HIGH)), list(chr, pos, peaknum, variant.ids)]

  magmean$sex <- ifelse(magmean$meanAltFreqLow > magmean$meanAltFreqHigh, "male", "female")


  seqSetFilter(genofile, variant.id=m$variant.ids)

  seqSetFilter(genofile, sample.id=scACids)

  msub <- m[, c("chr", "pos", "variant.ids", "sex", "rep", "deltaSNP"), with=TRUE]

    mrep1 <- msub[rep==1]
    mrep2 <- msub[rep==2]


# Pull out genotypes
  het <- t(seqGetData(genofile, "$dosage"))
  het <- as.data.table(het)

  colnames(het) <- c(seqGetData(genofile, "sample.id"))
  het$variant.ids <- seqGetData(genofile, "variant.id")

  setkey(het, variant.ids)
  setkey(snps, variant.ids)

  mhet <- merge(snps, het)

  mhetlong <- melt(mhet, measure.vars=scACids, variable.name="clone", value.name="dosage")

  setkey(scAC, clone)
  setkey(mhetlong, clone)
  mmhetlong <- merge(mhetlong, scAC)

  SCcounts <- mmhetlong[, .N, by=list(SC, variant.ids, dosage)]

#Remove NAs
  SCcounts <- SCcounts[dosage!="NA"]
  SCcountsA <- SCcounts[SC=="A"]
  SCcountsC <- SCcounts[SC=="C"]

  SCcountsAwide <- dcast(SCcountsA, variant.ids ~ dosage, value.var="N")
  colnames(SCcountsAwide) <- c("variant.ids", "Ados0", "Ados1", "Ados2")
  SCcountsAwide[is.na(Ados0),Ados0:=0]
  SCcountsAwide[is.na(Ados1),Ados1:=0]
  SCcountsAwide[is.na(Ados2),Ados2:=0]

  SCcountsCwide <- dcast(SCcountsC, variant.ids ~ dosage, value.var="N")
  colnames(SCcountsCwide) <- c("variant.ids", "Cdos0", "Cdos1", "Cdos2")
  SCcountsCwide[is.na(Cdos0),Cdos0:=0]
  SCcountsCwide[is.na(Cdos1),Cdos1:=0]
  SCcountsCwide[is.na(Cdos2),Cdos2:=0]

  setkey(SCcountsAwide, variant.ids)
  setkey(SCcountsCwide, variant.ids)
  SCcountsAC <- merge(SCcountsAwide, SCcountsCwide)

  SCcountsAC$Atot <- SCcountsAC$Ados0+SCcountsAC$Ados1+SCcountsAC$Ados2
  SCcountsAC$Ctot <- SCcountsAC$Cdos0+SCcountsAC$Cdos1+SCcountsAC$Cdos2


  setkey(SCcountsAC, variant.ids)
  setkey(mrep1, variant.ids)
  mmrep1 <- merge(SCcountsAC, mrep1)

  setkey(SCcountsAC, variant.ids)
  setkey(mrep2, variant.ids)
  mmrep2 <- merge(SCcountsAC, mrep2)

  setkey(SCcountsAC, variant.ids)
  setkey(magmean, variant.ids)
  mmagmean <- merge(SCcountsAC, magmean)


  mmrep1$AAlt <- round(mmrep1$Ados0*2+mmrep1$Ados1)/(mmrep1$Atot*2)
  mmrep1$CAlt <- round(mmrep1$Cdos0*2+mmrep1$Cdos1)/(mmrep1$Ctot*2)
  mmrep1$AaltminusCalt <- mmrep1$AAlt-mmrep1$CAlt

  mmrep2$AAlt <- round(mmrep2$Ados0*2+mmrep2$Ados1)/(mmrep2$Atot*2)
  mmrep2$CAlt <- round(mmrep2$Cdos0*2+mmrep2$Cdos1)/(mmrep2$Ctot*2)
  mmrep2$AaltminusCalt <- mmrep2$AAlt-mmrep2$CAlt

  mmagmean$AAlt <- round(mmagmean$Ados0*2+mmagmean$Ados1)/(mmagmean$Atot*2)
  mmagmean$CAlt <- round(mmagmean$Cdos0*2+mmagmean$Cdos1)/(mmagmean$Ctot*2)
  mmagmean$AaltminusCalt <- mmagmean$AAlt-mmagmean$CAlt

  mmrep1[,rAlt:=round(AAlt,1)]
  mmrep1[,rClt:=round(CAlt,1)]

  mmrep1[sex=="male", x:=((rClt) - (rAlt))]
  mmrep1[sex=="female", x:=((rAlt) - (rClt))]

  mmrep2[,rAlt:=round(AAlt,1)]
  mmrep2[,rClt:=round(CAlt,1)]

  mmrep2[sex=="male", x:=((rClt) - (rAlt))]
  mmrep2[sex=="female", x:=((rAlt) - (rClt))]

  mmagmean[,rAlt:=round(AAlt,1)]
  mmagmean[,rClt:=round(CAlt,1)]

  mmagmean[sex=="male", x:=((rClt) - (rAlt))]
  mmagmean[sex=="female", x:=((rAlt) - (rClt))]


  mmrep1$Cmale <- ifelse(mmrep1$sex=="male" & mmrep1$rClt==0.5, "1", ifelse(
    mmrep1$sex=="male" & mmrep1$rClt==1, "2", ifelse(
    mmrep1$sex=="female" & mmrep1$rClt==0.5, "-1", ifelse(
    mmrep1$sex=="female" & mmrep1$rClt==1, "-2", "0"
  ))))

  mmrep1$Amale <- ifelse(mmrep1$sex=="male" & mmrep1$rAlt==0.5, "1", ifelse(
    mmrep1$sex=="male" & mmrep1$rAlt==1, "2", ifelse(
    mmrep1$sex=="female" & mmrep1$rAlt==0.5, "-1", ifelse(
    mmrep1$sex=="female" & mmrep1$rAlt==1, "-2", "0"
  ))))



  setkey(mmrep1, chr, pos, variant.ids)
  setkey(snpsinpeaks, chr, pos, variant.ids)
  mmrep1B <- merge(mmrep1, snpsinpeaksnarrow)

  mmrep1B[x!=0,list(meanx=mean(x>0), N=length(x)), list(peaknum)]
  mmrep1B[,list(meanAmale=mean(as.numeric(Amale)), meanCmale=mean(as.numeric(Cmale))), list(peaknum)]


  setkey(mmrep2, chr, pos, variant.ids)
  setkey(snpsinpeaksnarrow, chr, pos, variant.ids)
  mmrep2B <- merge(mmrep2, snpsinpeaksnarrow)

  mmrep2B[x!=0,list(meanx=mean(x>0), N=length(x)), list(peaknum)]
  mmrep2B[,list(meanAmale=mean(as.numeric(Amale)), meanCmale=mean(as.numeric(Cmale))), list(peaknum)]

  setkey(mmagmean, chr, pos, variant.ids, peaknum)
  setkey(snpsinpeaksnarrow, chr, pos, variant.ids, peaknum)
  mmagmeanB <- merge(mmagmean, snpsinpeaksnarrow)

  mmagmeanB[x!=0,list(meanx=mean(x>0), N=length(x)), list(peaknum)]

  mmagmeanB_peak10 <- mmagmeanB[peaknum==10]
  peak10genes <- as.data.table(table(mmagmeanB_peak10$gene, mmagmeanB_peak10$col))
  colnames(peak10genes) <- c("gene", "type", "N")

  mmagmeanBgenes <- as.data.table(table(mmagmeanB$gene, mmagmeanB$col))
  colnames(mmagmeanBgenes) <- c("gene", "type", "N")

  ### gff
    gff <- fread("/scratch/kbb7sh/genomefiles/Daphnia.aed.0.6.gff")
    #gff[,ID:=gsub("ID=", "", first(tstrsplit(tstrsplit(V9, ";|:")[[1]], "-")))]
    gff[,ID:=gsub("ID=", "", tstrsplit(V9, ";|:")[[1]])]
    gff[,ss:=paste(V1, V3, V4, V5, ID, sep="_")]
    setkey(gff, ss)
    gff2 <- gff[!duplicated(gff)]
    dim(gff)
    dim(gff2)
    collapse_fun <- function(start, stop) {
      #start <- c(1,5); stop<-c(4, 10)
      length(unique(unlist(apply(cbind(start, stop), 1, function(x) x[1]:x[2]))))
    }
    gff.ag <- gff2[V3=="CDS",list(cds=sum(V5-V4 + 1), cds2=collapse_fun(start=V4, stop=V5)), list(chr=V1, id=ID)]
    gff.ag[,aa:=cds2/3]
    gff.ag
    summary(gff.ag)
    gff.ag[aa!=round(aa)]
    gff.ag[,gene:=first(tstrsplit(id, "-"))]
  ### use this
    gff.ag.ag <- gff.ag[,list(max_cds=max(cds), max_aa=max(aa)), list(gene)]

  ### Merge with peak10
    setkey(gff.ag.ag, gene)
    setkey(peak10genes, gene)
    mpeak10genes <- merge(peak10genes, gff.ag.ag)
    mpeak10genes$peraa <- mpeak10genes$N/mpeak10genes$max_aa
    ggplot(data=mpeak10genes[type=="missense_variant" | type=="synonymous_variant"], aes(x=peraa, fill=type)) + geom_histogram()

    mpeak10genespNpS <- mpeak10genes[type=="missense_variant" | type=="synonymous_variant"]
    mpeak10genespNpSwide <- dcast(mpeak10genespNpS, gene ~ type, value.var=c("N", "max_aa", "peraa"))
    mpeak10genespNpSwide$pN_pS <- mpeak10genespNpSwide$peraa_missense_variant/mpeak10genespNpSwide$peraa_synonymous_variant
    mpeak10genespNpSwide$totN <- mpeak10genespNpSwide$N_missense_variant+mpeak10genespNpSwide$N_synonymous_variant

    setkey(gff.ag.ag, gene)
    setkey(mmagmeanBgenes, gene)
    mmmagmeanBgenes <- merge(mmagmeanBgenes, gff.ag.ag)
    mmmagmeanBgenes$peraa <- mmmagmeanBgenes$N/mmmagmeanBgenes$max_aa

    mmmagmeanBgenespNpS <- mmmagmeanBgenes[type=="missense_variant" | type=="synonymous_variant"]
    mmmagmeanBgenespNpSwide <- dcast(mmmagmeanBgenespNpS, gene ~ type, value.var=c("N", "max_aa", "peraa"))
    mmmagmeanBgenespNpSwide$pN_pS <- mmmagmeanBgenespNpSwide$peraa_missense_variant/mmmagmeanBgenespNpSwide$peraa_synonymous_variant
    mmmagmeanBgenespNpSwide$totN <- mmmagmeanBgenespNpSwide$N_missense_variant+mmmagmeanBgenespNpSwide$N_synonymous_variant

    mmagmeanBsub <- mmagmeanB[, c("gene", "peaknum"), with=TRUE]
    mmagmeanBsubunique <- unique(mmagmeanBsub)
    setkey(mmagmeanBsubunique, gene)
    setkey(mmmagmeanBgenespNpSwide, gene)
    mmmagmeanBgenespNpSwidepeak <- merge(mmmagmeanBgenespNpSwide, mmagmeanBsubunique)

    save(mmmagmeanBgenespNpSwidepeak, file="mmmagmeanBgenespNpSwidepeak_20201124.Rdata")

  gff <- fread("/scratch/kbb7sh/genomefiles/Daphnia.aed.0.6.gff")
  colnames(gff) <- c("chr", "maker", "type", "start", "stop", "V6", "V7", "V8", "Notes")
  tmp <- sapply(strsplit(gff$Notes, ";"), getElement, 1)
  tmp2 <- sapply(strsplit(tmp, "="), getElement, 2)
  tmp3 <- sapply(strsplit(tmp2, ":"), getElement, 1)
  gff$gene <- tmp3

  gff$length <- gff$stop-gff$start


  uniprot <- fread("/scratch/kbb7sh/genomefiles/Daphnia_annotation_uniprot.txt")

  magmeanaddgenes <- foreach(snp.i=1:dim(mmagmeanB)[1], .combine="rbind")%do%{
    #snp.i <- 1
    c <- mmagmeanB$chr[[snp.i]]
    p <- mmagmeanB$pos[[snp.i]]

    tmp <- gff[chr==c & start <= p & stop >= p & type=="gene"]
    tmp2 <- mmagmeanB[chr==c & pos==p]
    tmp2$gene <-tmp$gene
    tmp2
  }


  allelefreqinpeaksnarrow <- foreach(peak.i=1:dim(peaks)[1], .combine="rbind")%do%{
    #peak.i <- 1
    c <- peaks$CHROM[[peak.i]]
    strt <- peaks$narrowstart[[peak.i]]
    stp <- peaks$narrowstop[[peak.i]]

    tmp <- alleleFreqs[chr==c & pos >= strt & pos <= stp]
    tmp$peaknum <- c(peak.i)
    tmp
  }

  setkey(allelefreqinpeaksnarrow, chr, pos, peaknum)
  setkey(snpsinpeaksnarrow, chr, pos, peaknum)
  m <- merge(snpsinpeaksnarrow, allelefreqinpeaksnarrow)

  snpsinpeaksnarrowstopgained <- snpsinpeaksnarrow[col=="stop_gained"]
  snpsinpeaksnarrowstopgainedids <- snpsinpeaksnarrowstopgained$variant.ids

  setkey(snpsinpeaksnarrowstopgained, chr, pos)
  setkey(alleleFreqs, chr, pos)
  msnpsinpeaksnarrowstopgained <- merge(snpsinpeaksnarrowstopgained, alleleFreqs)

  snpsinpeaksnarrowmissense_variant <- snpsinpeaksnarrow[col=="missense_variant"]
  snpsinpeaksnarrowmissense_variantids <- snpsinpeaksnarrowmissense_variant$variant.ids


  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  scAC <- sc[SC=="A" | SC=="C"]
  scACids <- scAC$clone
  seqSetFilter(genofile, sample.id=scACids)
  seqSetFilter(genofile, variant.id=snpsinpeaksnarrowmissense_variantids)

  # Pull out genotypes
  	het <- t(seqGetData(genofile, "$dose"))
  	het <- as.data.table(het)

  	colnames(het) <- c(seqGetData(genofile, "sample.id"))
  	het$variant.ids <- seqGetData(genofile, "variant.id")

  	setkey(het, variant.ids)
  	setkey(snps, variant.ids)

  	mhet <- merge(snps, het)

  	mhetlong <- melt(mhet, measure.vars=scACids, variable.name="clone", value.name="dosage")

  	setkey(scAC, clone)
  	setkey(mhetlong, clone)
  	mmhetlong <- merge(mhetlong, scAC)
  	mmhetlong <- mmhetlong[Nonindependent==0]

  	SCcounts <- mmhetlong[, .N, by=list(SC, variant.ids, dosage)]

    #Remove NAs
    	SCcounts <- SCcounts[dosage!="NA"]
    	SCcountsA <- SCcounts[SC=="A"]
    	SCcountsC <- SCcounts[SC=="C"]

    	SCcountsAwide <- dcast(SCcountsA, variant.ids ~ dosage, value.var="N")
    	colnames(SCcountsAwide) <- c("variant.ids", "Ados0", "Ados1", "Ados2")
    	SCcountsAwide[is.na(Ados0),Ados0:=0]
    	SCcountsAwide[is.na(Ados1),Ados1:=0]
    	SCcountsAwide[is.na(Ados2),Ados2:=0]

    	SCcountsCwide <- dcast(SCcountsC, variant.ids ~ dosage, value.var="N")
    	colnames(SCcountsCwide) <- c("variant.ids", "Cdos0", "Cdos1", "Cdos2")
    	SCcountsCwide[is.na(Cdos0),Cdos0:=0]
    	SCcountsCwide[is.na(Cdos1),Cdos1:=0]
    	SCcountsCwide[is.na(Cdos2),Cdos2:=0]

    	setkey(SCcountsAwide, variant.ids)
    	setkey(SCcountsCwide, variant.ids)
    	SCcountsAC <- merge(SCcountsAwide, SCcountsCwide)

      SCcountsAC$totA <- SCcountsAC$Ados0+SCcountsAC$Ados1+SCcountsAC$Ados2
      SCcountsAC$totC <- SCcountsAC$Cdos0+SCcountsAC$Cdos1+SCcountsAC$Cdos2
      SCcountsAC$prophetA <- SCcountsAC$Ados1/SCcountsAC$totA
      SCcountsAC$prophetC <- SCcountsAC$Cdos1/SCcountsAC$totC
      SCcountsAC$hetAminushetC <- abs(SCcountsAC$prophetA-SCcountsAC$prophetC)




### Stop gained mutations in narrow peaks
variant.ids Ados0 Ados1 Ados2 Cdos0 Cdos1 Cdos2
1:      315442     0    85     0     0    28     0
2:      610456     0    63     2     0    17     2
3:      610468     0     2    86     0     0    28
4:      610504     0     0    90     0     0    28
5:      610509    88     1     0    28     0     0
6:      610518     0    87     1     0    28     0
7:      610520     0    87     1     0    28     0
8:      611040     0     0    86     0     0    28
