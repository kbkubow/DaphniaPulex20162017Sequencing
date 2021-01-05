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

### Load in gff file
  gff <- fread("/scratch/kbb7sh/genomefiles/Daphnia.aed.0.6.gff")
  colnames(gff) <- c("chr", "maker", "type", "start", "stop", "V6", "V7", "V8", "Notes")
  tmp <- sapply(strsplit(gff$Notes, ";"), getElement, 1)
  tmp2 <- sapply(strsplit(tmp, "="), getElement, 2)
  tmp3 <- sapply(strsplit(tmp2, ":"), getElement, 1)
  gff$gene <- tmp3

  uniprot <- fread("/scratch/kbb7sh/genomefiles/Daphnia_annotation_uniprot.txt")

### Lets look at the peak on 2217 at pos 3,968,943. Only have top SNP right now, lets look in 5000 window up and downstream
  chr2217peak <- gff[chr=="Scaffold_2217_HRSCAF_2652" & start > 3963943 & stop < 3973943]
  chr2217peakgenes <- data.table(qseqid=unique(chr2217peak$gene))
  setkey(chr2217peakgenes, qseqid)
  setkey(uniprot, qseqid)
  mchr2217peakgenes <- merge(chr2217peakgenes, uniprot)

  # Two genes appear in this stretch, mostly Daphnia00600-RD, some Daphnia00600-RB, Daphnia00600-RA,
  # Daphnia00600-RE, Daphnia00600-RC, and small bit of Daphnia00601-RA
  # Daphnia00600-RD -> tr|E9GCB6|E9GCB6_DAPPU - no panther protein class or family subfamily
  # Daphnia00600-RA -> tr|E9GCB6|E9GCB6_DAPPU
  # Daphnia00600-RB -> tr|E9GCB7|E9GCB7_DAPPU - Metalloendopeptidase - PTHR45645:SF23
  # Daphnia00600-RE -> tr|E9GCB7|E9GCB7_DAPPU
  # Daphnia00600-RC -> tr|E9GCB7|E9GCB7_DAPPU
  # Daphnia00601-RA -> tr|E9GCB4|E9GCB4_DAPPU - Lipase - PTHR11005:FS100

  chr7757peak <- gff[chr=="Scaffold_7757_HRSCAF_8726" & start > 7873074 & stop < 9349506]
  chr7757peakgenes <- data.table(qseqid=unique(chr7757peak$gene))
  setkey(chr7757peakgenes, qseqid)
  setkey(uniprot, qseqid)
  mchr7757peakgenes <- merge(chr7757peakgenes, uniprot)

  chr9199peak <- gff[chr=="Scaffold_9199_HRSCAF_10755" & start > 4964654 & stop < 6547651]
  chr9199peakgenes <- data.table(qseqid=unique(chr9199peak$gene))
  setkey(chr9199peakgenes, qseqid)
  setkey(uniprot, qseqid)
  mchr9199peakgenes <- merge(chr9199peakgenes, uniprot)



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
  # Collapsing additional annotations to original SNP vector length
      snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=",")),
                        list(variant.ids=id)]
      snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]

### Add annotations to snp table
  setkey(snps, variant.ids)
  setkey(snp.dt1.an, variant.ids)
  snpsann <- merge(snps, snp.dt1.an)

  # Load in peaks file
    #peaks <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/peaks.csv")
    load("gprime_peaks.replicates.250K.05.Rdata")

    peaks$narrowstart <- peaks$posMaxGprime-20000
    peaks$narrowstop <- peaks$posMaxGprime+20000


  # Figure out which snps fall in peaks

  snpsinpeaks <- foreach(peak.i=1:dim(peaks)[1], .combine="rbind")%do%{
    #peak.i <- 1
    c <- peaks$chr[[peak.i]]
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

    tmp <- snpsann[chr==c & pos >= strt-125000 & pos <= stp+125000]
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

  alleleFreqs$ALT_FREQ_LOW <- alleleFreqs$AD_ALT.LOW/alleleFreqs$DP.LOW
  alleleFreqs$ALT_FREQ_HIGH <- alleleFreqs$AD_ALT.HIGH/alleleFreqs$DP.HIGH
  alleleFreqs$ALT_FREQ_HIGH_MINUS_LOW <- (alleleFreqs$ALT_FREQ_HIGH-alleleFreqs$ALT_FREQ_LOW)

  setkey(snpsinpeaksnarrow, chr, pos)
  setkey(alleleFreqs, chr, pos)
  m <- merge(snpsinpeaksnarrow, alleleFreqs)

  m$sex <- ifelse(m$ALT_FREQ_LOW > m$ALT_FREQ_HIGH, "male", "female")

  mrep1 <- m[rep==1]



  seqSetFilter(genofile, variant.id=m$variant.ids)

  seqSetFilter(genofile, sample.id=scACids)

  msub <- m[, c("chr", "pos", "variant.ids", "sex"), with=TRUE]

  mag <- m[,list(avgALT_FREQ_LOW = mean(ALT_FREQ_LOW, na.rm=TRUE), avgALT_FREQ_HIGH=
    mean(ALT_FREQ_HIGH, na.rm=TRUE)), list(chr,pos,variant.ids) ]

  mag$sex <- ifelse(mag$avgALT_FREQ_LOW>mag$avgALT_FREQ_HIGH, "male", "female")


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


  msubunique <- unique(msub)

  setkey(SCcountsAC, variant.ids)
  setkey(mag, variant.ids)
  m2 <- merge(SCcountsAC, mag)

  m2$Asex <- ifelse(m2$sex=="male" & m2$Ados1 >= m2$Atot-3, "male", ifelse(m2$sex=="female" &
    m2$Ados1 >=m2$Atot-3, "female", ifelse(m2$sex=="male" & m2$Ados1 < m2$Atot-3, "female", "male")))
  m2$Csex <- ifelse(m2$sex=="male" & m2$Cdos1 >= m2$Ctot-3, "male", ifelse(m2$sex=="female" &
    m2$Cdos1 >=m2$Ctot-3, "female", ifelse(m2$sex=="male" & m2$Cdos1 < m2$Ctot-3, "female", "male")))

  m2$propAhet <- m2$Ados1/m2$Atot
  m2$propChet <- m2$Cdos1/m2$Ctot
  m2$absdiffhet <- abs(m2$propAhet-m2$propChet)

  setkey(m2, chr, pos, variant.ids)
  setkey(snpsinpeaksnarrow, chr, pos, variant.ids)
  m3 <- merge(m2, snpsinpeaksnarrow)


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
