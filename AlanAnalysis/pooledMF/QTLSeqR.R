### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R


#install.packages("devtools")
#
## use devtools to install QTLseqr
#devtools::install_github("bmansfeld/QTLseqr")


### libraries
  library(data.table)
  library(QTLseqr)

#### using the ASE read counter data
### load data
  male <- fread("/mnt/sammas_storage/bergland-lab/alan/DBunkMale.pooledAF.aseReadCounter.allvariant.delim")
  pe <- fread("/mnt/sammas_storage/bergland-lab/alan/DBunkPE.pooledAF.aseReadCounter.allvariant.delim")

  setnames(male,
           c("contig", "position", "refAllele", "altAllele", "refCount", "altCount"),
           c("CHROM", "POS", "REF", "ALT", "AD_REF.LOW", "AD_ALT.LOW"))


   setnames(pe,
            c("contig", "position", "refAllele", "altAllele", "refCount", "altCount"),
            c("CHROM", "POS", "REF", "ALT", "AD_REF.HIGH", "AD_ALT.HIGH"))


  setkey(male, CHROM, POS)
  setkey(pe, CHROM, POS)

  m <- merge(male[,c("CHROM", "POS", "REF", "ALT", "AD_REF.LOW", "AD_ALT.LOW"),with=F],
             pe[,c("CHROM", "POS", "REF", "ALT", "AD_REF.HIGH", "AD_ALT.HIGH"),with=F])

  m[,deltaSNP:=qlogis(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - qlogis(AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
  m.ag <- m[,.N,CHROM]
  m <-  m[CHROM%in%m.ag[N>1000]$CHROM]

  m[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
  m[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
  m[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]

  ### filter?
  df_filt <- filterSNPs(SNPset = as.data.frame(m),
                       refAlleleFreq = 0.1,
                       minTotalDepth = 100,
                       maxTotalDepth = 1500,
                       depthDifference = 200,
                       minSampleDepth = 100,
                       verbose = TRUE)


  gprime <- runGprimeAnalysis(SNPset=df_filt, windowSize=200000)
  gprime <- as.data.table(gprime)
  p <- plotQTLStats(SNPset = gprime, var = "Gprime") +
       geom_hline(aes(yintercept=min(gprime[qvalue<.1]$Gprime)), color="red") +
       geom_hline(aes(yintercept=min(gprime[qvalue<.05]$Gprime)), color="blue") +
       geom_hline(aes(yintercept=min(gprime[qvalue<.01]$Gprime)), color="green")
  p
  peaks <- as.data.table(getQTLTable(SNPset = gprime, method = "Gprime", alpha = 0.01, export = FALSE))
  peaks[order(meanGprime)]



### VarScan data
### libraries
  library(data.table)
  library(QTLseqr)
  library(ggplot2)

### pooled data
  load("/mnt/sammas_storage/bergland-lab/alan/totalADRDlongall.Rdata")

  geno[,effRD:=floor(RRD)]
  geno[,effPA:=round(propalt*effRD)/effRD]


  #### DBunk
    geno.w <- dcast(geno[pond=="DBunk"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

    geno.w.ag <- geno.w[,.N,chr]
    geno.w <- geno.w[chr%in%geno.w.ag[N>1000]$chr]

    geno.w.comb <- geno.w[,list(AD_REF.LOW = (1-effPA_DBunkMale1)*effRD_DBunkMale1 + (1-effPA_DBunkMale2)*effRD_DBunkMale2,
                                AD_ALT.LOW = (effPA_DBunkMale1)*effRD_DBunkMale1 + (effPA_DBunkMale2)*effRD_DBunkMale2,
                                AD_REF.HIGH = (1-effPA_DBunkPE1)*effRD_DBunkPE1 + (1-effPA_DBunkPE2)*effRD_DBunkPE2,
                                AD_ALT.HIGH = (effPA_DBunkPE1)*effRD_DBunkPE1 + (effPA_DBunkPE2)*effRD_DBunkPE2),
                            list(CHROM=chr, POS=pos)]

    geno.w.comb[,deltaSNP:=(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - (AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
    geno.w.comb[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
    geno.w.comb[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
    geno.w.comb[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]

    geno.w.comb.filt <- filterSNPs(SNPset = as.data.frame(geno.w.comb),
                         refAlleleFreq = 0.05,
                         minTotalDepth = 50,
                         maxTotalDepth = 350,
                         depthDifference = 200,
                         minSampleDepth = 25,
                         verbose = TRUE)


    gprime <- runGprimeAnalysis(SNPset=geno.w.comb.filt, windowSize=250000, outlierFilter="Hampel", maxk=100)
    gprime <- as.data.table(gprime)
    p <- plotQTLStats(SNPset = gprime, var = "Gprime") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.1]$Gprime)), color="red") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.05]$Gprime)), color="blue") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.01]$Gprime)), color="green") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.001]$Gprime)), color="green")

    p
    peaks <- as.data.table(getQTLTable(SNPset = gprime, method = "Gprime", alpha = 0.0001, export = FALSE))
    peaks[order(meanGprime)]


  ### DBunk
    geno.w <- dcast(geno[pond=="DBunk"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

    geno.w.ag <- geno.w[,.N,chr]
    geno.w <- geno.w[chr%in%geno.w.ag[N>1000]$chr]

    geno.w.comb <- geno.w[,list(AD_REF.LOW = (1-effPA_DBunkMale)*effRD_DBunkMale,
                                AD_ALT.LOW = (effPA_DBunkMale)*effRD_DBunkMale,
                                AD_REF.HIGH = (1-effPA_DBunkPE1)*effRD_DBunkPE1 + (1-effPA_DBunkPE2)*effRD_DBunkPE2,
                                AD_ALT.HIGH = (effPA_DBunkPE1)*effRD_DBunkPE1 + (effPA_DBunkPE2)*effRD_DBunkPE2),
                            list(CHROM=chr, POS=pos)]

    geno.w.comb[,deltaSNP:=(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - (AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
    geno.w.comb[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
    geno.w.comb[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
    geno.w.comb[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]



    geno.w.comb.filt <- filterSNPs(SNPset = as.data.frame(geno.w.comb),
                         refAlleleFreq = 0.05,
                         minTotalDepth = 50,
                         maxTotalDepth = 350,
                         depthDifference = 200,
                         minSampleDepth = 25,
                         verbose = TRUE)


    gprime <- runGprimeAnalysis(SNPset=geno.w.comb.filt, windowSize=250000, outlierFilter="Hampel", maxk=100)
    gprime <- as.data.table(gprime)
    p <- plotQTLStats(SNPset = gprime, var = "Gprime") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.1]$Gprime)), color="red") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.05]$Gprime)), color="blue") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.01]$Gprime)), color="green") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.001]$Gprime)), color="green")

    p
    peaks <- as.data.table(getQTLTable(SNPset = gprime, method = "Gprime", alpha = 0.0001, export = FALSE))
    peaks[order(meanGprime)]

  ### DCat
    geno.w <- dcast(geno[pond=="DCat"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

    geno.w.ag <- geno.w[,.N,chr]
    geno.w <- geno.w[chr%in%geno.w.ag[N>1000]$chr]

    geno.w.comb <- geno.w[,list(AD_REF.LOW = (1-effPA_DCatMale)*effRD_DCatMale,
                                AD_ALT.LOW = (effPA_DCatMale)*effRD_DCatMale,
                                AD_REF.HIGH = (1-effPA_DCatPE1)*effRD_DCatPE1 + (1-effPA_DCatPE2)*effRD_DCatPE2,
                                AD_ALT.HIGH = (effPA_DCatPE1)*effRD_DCatPE1 + (effPA_DCatPE2)*effRD_DCatPE2),
                            list(CHROM=chr, POS=pos)]

    geno.w.comb[,deltaSNP:=(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - (AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
    geno.w.comb[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
    geno.w.comb[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
    geno.w.comb[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]



    geno.w.comb.filt <- filterSNPs(SNPset = as.data.frame(geno.w.comb),
                         refAlleleFreq = 0.05,
                         minTotalDepth = 50,
                         maxTotalDepth = 350,
                         depthDifference = 200,
                         minSampleDepth = 25,
                         verbose = TRUE)



    gprime <- runGprimeAnalysis(SNPset=geno.w.comb.filt, windowSize=250000, outlierFilter="Hampel", maxk=100)
    gprime <- as.data.table(gprime)
    p <- plotQTLStats(SNPset = gprime, var = "Gprime") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.1]$Gprime)), color="red") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.05]$Gprime)), color="blue") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.01]$Gprime)), color="green") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.001]$Gprime)), color="green")

    p
    peaks <- as.data.table(getQTLTable(SNPset = gprime, method = "Gprime", alpha = 0.0001, export = FALSE))
    peaks[order(meanGprime)]









ggplot(data=gprime, aes(x=POS, y=Gprime)) + geom_line() + facet_wrap(~CHROM, nrow=1)

getQTLTable(SNPset = gprime, method = "Gprime",alpha = 0.005, export = FALSE)


gprime <- as.data.table(gprime)
gprime[,ed4:=sqrt(deltaSNP^2 + (-1*deltaSNP)^2)^4]
ggplot(data=gprime, aes(x=POS, y=ed4)) + geom_line() + facet_wrap(~CHROM, nrow=1)

o <- runQTLseqAnalysis(df_filt, windowSize = 500000,
popStruc = "F2",
bulkSize = c(70, 100), replications = 10000,
intervals = c(95, 99) )

plotQTLStats(SNPset = o, var = "deltaSNP", plotIntervals = TRUE))


 p3 <- plotQTLStats(SNPset = gprime, var = "Gprime", plotThreshold = TRUE, q = 0.1)
 p3


 gprime <- as.data.table(gprime)
 gprime[qvalue<.01]
