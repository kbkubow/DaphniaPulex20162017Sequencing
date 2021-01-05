
### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(4)

#### using the ASE read counter data
### load data
  male <- fread("~/D8Male.pooledAF.aseReadCounter.allvariant.delim")
  pe <- fread("~/D8PE.pooledAF.aseReadCounter.allvariant.delim")

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

  m[,male.freq:=AD_REF.LOW/DP.LOW]
  m[,pe.freq:=AD_REF.HIGH/DP.HIGH]

### load in phased
  trio <- seqOpen("/Users/alanbergland/Desktop/testTrio.phased.gds")
  trio.dt <- data.table(CHROM=seqGetData(trio, "chromosome"),
                        POS=seqGetData(trio, "position"))
  seqSetFilter(trio, sample.id="April_2017_D8_103")
  trio.dt[,A.geno:=seqAlleleFreq(trio)]
  seqSetFilter(trio, sample.id="April_2017_D8_125")
  trio.dt[,B.geno:=seqAlleleFreq(trio)]

### combine
  setkey(m, CHROM, POS)
  setkey(trio.dt, CHROM, POS)
  mm <- merge(m, trio.dt)
  mm[,exp.freq:=A.geno/2 + B.geno/2]
  mm.ag <- na.omit(mm[,list(mu=c(mean(male.freq), mean(pe.freq)),
                            sd=c(sd(male.freq), sd(pe.freq)),
                            pool=c("male", "pe")),
              list(A.geno, B.geno)])




### generate windowed alelle frequency estimates

  chrs <- mm[, list(min=min(POS), max=max(POS), .N), CHROM]
  window.bp <- 200000
  step.bp <- 5000 # 5000

  wins <- foreach(i=1:dim(chrs)[1], .combine="rbind")%dopar%{
    #i<-1
    print(i)
    tmp <- data.table(chr=chrs[i]$CHROM, start=seq(from=chrs[i]$min, to=chrs[i]$max-window.bp, by=step.bp))
    tmp[,stop:=start+window.bp]
    tmp
  }
  wins[,i:=1:dim(wins)[1]]


  i <- wins[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(start-8660157))]$i

  maf <- 0.05
  o <- foreach(i=c(1:dim(wins)[1]), .errorhandling="remove")%do%{
    if(i%%100==0) print(paste(i, dim(wins)[1]))
    i <- 20377
    tmp <- mm[CHROM==wins$chr[i] & POS>=wins$start[i] & POS<=wins$stop[i]][male.freq>=maf & male.freq<=(1-maf) & pe.freq>=maf & pe.freq<=(1-maf)]
    ggplot(na.omit(tmp), aes(x=male.freq, y=pe.freq, color=interaction(A.geno, B.geno))) + geom_point() + ylim(0,1) + xlim(0,1) + geom_abline(aes(intercept=0,slope=1))


    tmp.ag <- (na.omit(tmp[,list(delta.mu=mean(deltaSNP), delta.sd=sd(deltaSNP), i=i, .N), list(A.geno, B.geno)]))
    return(tmp.ag)
  }

  o <- rbindlist(o)
  o.orig <- o

  o <- o[!(A.geno==1 & B.geno==1)]
  o <- o[!(A.geno==0 & B.geno==0)]

  o[,p.upper:=pnorm(0, delta.mu, delta.sd)]
  o[,p.lower:=1-pnorm(0, delta.mu, delta.sd)]
  o[,p:=2*apply(cbind(o$p.upper, o$p.lower), 1, min)]

  o.ag <- o[N>20,list(wm=weighted.mean(abs(delta.mu), w=N), N=sum(N)), list(i)]



  ggplot(o.ag, aes(x=i, y=wm)) + geom_line() + geom_point()
  o[i==20379]

  o.ag[wm>.68]












  o.ag <- o[,list(t=(mu.x[pool=="male"] - mu.x[pool=="pe"])/sqrt(sd.x[pool=="male"]^2 + sd.x[pool=="pe"]^2), n=mean(N)),
     list(A.geno, B.geno, i)]
  o.ag[,p:=1-pt(t, n-2)]

  summary(o.ag[n>50]$t)


  ggplot(data=o, aes(x=i, y=-log10(p))) + geom_line() +


  ggplot(data=tmp, aes(x=male.freq, y=pe.freq)) + geom_point() + ylim(0,1) + xlim(0,1)



  abline(0,1)


  lo <- loess(male.freq~POS, m)


  hist(m[male.freq<.95 & male.freq>.05]$male.freq, breaks=1000)





getG <- function(LowRef, HighRef, LowAlt, HighAlt)
{
    exp <- c(
        (LowRef + HighRef) * (LowRef + LowAlt) / (LowRef + HighRef + LowAlt + HighAlt),
        (LowRef + HighRef) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
        (LowRef + LowAlt) * (LowAlt + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
        (LowAlt + HighAlt) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt)
    )
    obs &lt;- c(LowRef, HighRef, LowAlt, HighAlt)

    G <
        2 * (rowSums(obs * log(
            matrix(obs, ncol = 4) / matrix(exp, ncol = 4)
        )))
    return(G)
}
