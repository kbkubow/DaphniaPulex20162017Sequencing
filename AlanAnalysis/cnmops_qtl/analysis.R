library(data.table)
library(foreach)

### load sample data
  samps <- fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### load cnmops
  load("/Users/alanbergland/cnmopsregionslong.Rdata")

  dat <- as.data.table(cnmopsregionslong)
  dat[,N:=as.numeric(gsub("CN", "", ploidy))]

  dat[,locus:=paste(seqnames, start, end)]

  dat <- merge(dat, samps, by="clone")

  dat.ag <- dat[SC%in%c("A", "C"), list(A.N= median(N[SC=="A"]), C.N=median(N[SC=="C"])), list(locus)]

  table(dat$locus)

  dat.ag[]


### load phenotype data
  load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")

### for each locus
  loci <- dat.ag[(A.N!=C.N) | (A.N==1 & C.N==1)]$locus

  o <- foreach(i=loci, .errorhandling="remove")%do%{
    # i<-loci[1]
    #i <- "Scaffold_2217_HRSCAF_2652 4405501 4427250"
    message(i)
    tmp <- dat[locus==i]
    tmp <- merge(tmp, male, by="clone")

    tmp.ag <- tmp[,list(nUniq=length(unique(N)), N=median(N), males=sum(Males), tot=sum(NewTotal)), list(SCB, gr)]

    t1 <- glm(I(males/tot)~N, tmp.ag[gr=="AxC"], weights=tot, family=binomial())
    #t1 <- lm(propmale~N, tmp[gr=="AxC"])

    tmp.ag[gr=="AxC"][,pred:=predict(t1)]
  #  ggplot(data=tmp.ag[gr=="AxC"], aes(x=N, y=qlogis(males/tot))) + geom_point() + geom_abline(slope=summary(t1)$coef[2,1], intercept=summary(t1)$coef[1,1])
  ggplot(data=tmp.ag[gr=="AxC"], aes(x=N, y=(males/tot))) +
  geom_point() +
  geom_abline(slope=summary(t1)$coef[2,1], intercept=summary(t1)$coef[1,1])


    data.table(locus=i, p=summary(t1)$coef[2,4], beta=summary(t1)$coef[2,1], A=dat.ag[locus==i]$A.N, C=dat.ag[locus==i]$C.N)

  }
  o <- rbindlist(o)

  o[,chr:=tstrsplit(locus, " ")[[1]]]
  o[,start:=tstrsplit(locus, " ")[[2]]]
  o[,end:=tstrsplit(locus, " ")[[3]]]
  o[,mid:=as.numeric(end)/2 + as.numeric(start)/2]
  o[,concord:=sign(beta)==sign()]

  fisher.test(table(sign(o[p<.005]$beta), sign(o[p<.005]$A-o[p<.005]$C)))
  ggplot(data=o, aes(x=mid, y=-log10(p), color=chr)) + geom_line() + facet_grid(~chr)


  o[which.min(p)]
