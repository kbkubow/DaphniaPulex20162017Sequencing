### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)
  library(SeqArray)

### load genotyping data
  #load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
  load(file="/scratch/aob2x/daphnia_hwe_sims/subFiles.Rdata")
  sc[,year:=tstrsplit(clone, "_")[[2]]]

  ### open GDS object
  #genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds", allow.duplicate=TRUE)
  #genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
  genofile <- seqOpen("/scratch/aob2x/daphnia_hwe_sims/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### get differences between A&B
  genodat <- foreach(sc.i=c("A", "B"))%do%{
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sc[SC==sc.i][pond=="D8"]$clone)

    tmp <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     variant.id=seqGetData(genofile, "variant.id"),

                     sc=sc.i,
                     af= 1 - seqAlleleFreq(genofile))
  }
  genodat <- rbindlist(genodat)

  genodat.w <- dcast(genodat, chr + pos + variant.id ~ sc)


### pooled data
  load("/nv/vol186/bergland-lab/alan/totalADRDlongall.Rdata")

  geno[,effRD:=floor(RRD)]
  geno[,effPA:=round(propalt*effRD)/effRD]

  geno.w <- dcast(geno[pond=="D8"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

### merge
  setkey(genodat.w, chr, pos)
  setkey(geno.w, chr, pos)

  m <- merge(genodat.w, geno.w)

  m[,BA_delta:=B-A]
  m[,mf_delta.naive:=(effPA_D8Male1+effPA_D8Male2)/2 - (effPA_D8PE1+effPA_D8PE2)/2]
  m[,m.hat:=(effPA_D8Male1 * effRD_D8Male1 + effPA_D8Male2 * effRD_D8Male2) / (effRD_D8Male1 + effRD_D8Male2)]
  m[,f.hat:=(effPA_D8PE1 * effRD_D8PE1 + effPA_D8PE2 * effRD_D8PE2) / (effRD_D8PE1 + effRD_D8PE2)]
  m[,mf_delta:=m.hat - f.hat]


  setkey(snp.dt, chr, pos)
  setkey(m, chr, pos)
  m <- merge(snp.dt, m)


  m[]
  save(m, file="/nv/vol186/bergland-lab/alan/pool_AB.Rdata")

  table(sign(m[BA_delta!=0 & mf_delta!=0 & final.use==T]$BA_delta), sign(m[BA_delta!=0 & mf_delta!=0 & final.use==T]$mf_delta))

  fisher.test(table(sign(m[BA_delta!=0 & abs(mf_delta)>.25 & final.use==T]$BA_delta), sign(m[BA_delta!=0 & abs(mf_delta)>.25 & final.use==T]$mf_delta)))
  fisher.test(table(sign(m[BA_delta!=0 & abs(mf_delta)>.25 & use.chr==T]$BA_delta), sign(m[BA_delta!=0 & abs(mf_delta)>.25 & use.chr==T]$mf_delta)))
  fisher.test(table(sign(m[BA_delta!=0 & abs(mf_delta)>.25]$BA_delta), sign(m[BA_delta!=0 & abs(mf_delta)>.25]$mf_delta)))

  fisher.test(table(sign(m[BA_delta!=0 & abs(mf_delta)>.25 & use.chr==T]$BA_delta), sign(m[BA_delta!=0 & abs(mf_delta)>.25 & use.chr==T]$mf_delta)))



  m.ag <- foreach(d=seq(from=0, to=.25, by=.01), .combine="rbind")%do%{
    print(d)
    #d <- .05
    fet <- fisher.test(table(sign(m[BA_delta!=0 & abs(mf_delta)>d & use.chr==T]$BA_delta),
                  sign(m[BA_delta!=0 & abs(mf_delta)>d & use.chr==T]$mf_delta)))

    data.table(delta=d, or=fet$estimate, p=fet$p.value, lci=fet$conf.int[1], uci=fet$conf.int[2])
  }

  d <- .25
  m[BA_delta!=0 & abs(mf_delta)>d & use.chr==T][,c("chr", "A", "B", "m.hat", "f.hat"), with=F]

  save(m, m.ag, file="/nv/vol186/bergland-lab/alan/pool_AB.Rdata")







###################################
#### here

### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(doMC)
  registerDoMC(20)

### load data from Rivanna
  load("/mnt/sammas_storage/bergland-lab/alan/pool_AB.Rdata")
  m[!is.na(A),A.geno := unlist(sapply(m[!is.na(A)]$A, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[!is.na(B),B.geno := unlist(sapply(m[!is.na(B)]$B, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[,BA.geno.diff := B.geno - A.geno]
  m[,inform:=!((A.geno==0 & B.geno==0) | (A.geno==2 & B.geno==2))]


### generate weighted allele frequencies from pools
  m[,D8Male.f:=(effRD_D8Male1*effPA_D8Male1 + effRD_D8Male2*effPA_D8Male2) / (effRD_D8Male1+effRD_D8Male2)]
  m[,D8Male.rd:=(effRD_D8Male1+effRD_D8Male2)]

  m[,D8PE.f:=(effRD_D8PE1*effPA_D8PE1 + effRD_D8PE2*effPA_D8PE2) / (effRD_D8PE1+effRD_D8PE2)]
  m[,D8PE.rd:=(effRD_D8PE1+effRD_D8PE2)]

### polarize
  m[,A.geno.polarize := A.geno]
  m[,B.geno.polarize := B.geno]

  m[A.geno==2, A.geno.polarize := 0]
  m[A.geno==2 & B.geno==0, B.geno.polarize := 2]
  m[A.geno==2 & B.geno==2, B.geno.polarize := 0]
  m[,BA.geno.diff.polarize := B.geno.polarize - A.geno.polarize]

  m[,D8Male.f.polarize:=D8Male.f]
  m[,D8PE.f.polarize:=D8PE.f]
  m[A.geno==2, D8Male.f.polarize := 1 - D8Male.f]
  m[A.geno==2, D8PE.f.polarize := 1 - D8PE.f]

### fold
  m[,A.geno.fold := A.geno]
  m[,B.geno.fold := B.geno]
  m[A.geno==2, A.geno.fold := 0]
  m[B.geno==2, B.geno.fold := 0]

  m[,D8Male.f.fold:=D8Male.f]
  m[,D8PE.f.fold:=D8PE.f]
  m[D8Male.f>.5, D8Male.f.fold:= 1- D8Male.f]
  m[D8PE.f>.5, D8PE.f.fold := 1 - D8PE.f]

  m[,BA.geno.diff.fold := B.geno.fold - A.geno.fold]


### simple odds ratio
  m[, or.polarize:=((D8Male.f.polarize*D8Male.rd)/((1-D8Male.f.polarize)*D8Male.rd)) / ((D8PE.f.polarize*D8PE.rd)/((1-D8PE.f.polarize)*D8PE.rd))]
  m[, or:=((D8Male.f*D8Male.rd)/((1-D8Male.f)*D8Male.rd)) / ((D8PE.f*D8PE.rd)/((1-D8PE.f)*D8PE.rd))]
  m[, or.fold:=((D8Male.f.fold*D8Male.rd)/((1-D8Male.f.fold)*D8Male.rd)) / ((D8PE.f.fold*D8PE.rd)/((1-D8PE.f.fold)*D8PE.rd))]

### filter down a small bit
  m.inform <- m[final.use==T][inform==T][or!=0 & or!=Inf]
  setkey(m.inform, chr, pos)

### generate windows
  chrs <- m.inform[use.chr==T, list(min=min(pos), max=max(pos)), chr]
  window.bp <- 50000
  step.bp <- 5000

  wins <- foreach(i=1:dim(chrs)[1], .combine="rbind")%dopar%{
    #i<-1
    print(i)
    tmp <- data.table(chr=chrs[i]$chr, start=seq(from=chrs[i]$min, to=chrs[i]$max-window.bp, by=step.bp))
    tmp[,stop:=start+window.bp]
    tmp
  }
  wins[,i:=1:dim(wins)[1]]

  setkey(m.inform, chr, pos)

### run SW analysis

  o <- foreach(i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
    #i<-145167 ; i <- 58987
    #i <- wins[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(start-8660157))]$i
    #i <- wins[chr=="Scaffold_6786_HRSCAF_7541"][which.min(abs(start-12960054))]$i
    if(i%%100==0) print(paste(i, dim(wins)[1], sep=" / "))


    m.inform[,win:=0]
    m.inform[chr==wins$chr[i] & pos>=wins$start[i] & pos<=wins$stop[i], win:=1]

    mod <- summary(lm(log2(or)~((BA.geno.diff))*win, m.inform))

    data.table(win.i=i, chr=wins$chr[i], start=wins$start[i], stop=wins$stop[i], nSNPs=sum(m.inform$win),
              r2=mod$r.squared,
              beta.int=mod$coefficients[4,1],
              se=mod$coefficients[4,2],
              p=(mod$coefficients[4,4]))


  #data.table(win.i=i, chr=wins$chr[i], start=wins$start[i], stop=wins$stop[i], nSNPs=sum(m.inform$win),
  #         p=anova(mod)$P[3])

    }

  ### save
    o <- rbindlist(o)
    o[,pa:=p.adjust(p, "bonferroni")]
    save(o, file="~/slidingWindow_oddsRatio.Rdata")

  ### load
    load(file="~/slidingWindow_oddsRatio.Rdata")

  ### polarized Manhattan plot
    o[,ppa:=(sign(beta.int)*log10(pa))]

    gwas.plot <- ggplot(data=o[nSNPs>25], aes(y=(ppa), x=win.i, color=chr)) + geom_line()

  ### what do these windows look like
    plot.it <- function(i=o[which.min(ppa)]$win.i) {
      #i <- wins[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(start-8660157))]$i
      #i <- o[nSNPs>25][which.max(ppa)]$win.i
      #i <- o[win.i<5000][which.min(ppa)]$win.i

      m.inform[,win:=0]
      m.inform[chr==wins$chr[i] & pos>=wins$start[i] & pos<=wins$stop[i], win:=1]

      ml <- melt(m.inform[,c("chr", "pos", "win", "D8Male.f", "D8PE.f", "BA.geno.diff", "A.geno", "B.geno"), with=F],
                 id.vars=c("chr", "pos", "win", "BA.geno.diff", "A.geno", "B.geno"))
      ml[,winf:=c("genome", "targetWindow")[win+1]]
      ml[,exp_hom:= 1 - 2*value*(1-value)]

      or.plot <- ggplot(m.inform, aes(x=as.factor((BA.geno.diff)), y=log2(or.fold),
                           group=interaction((BA.geno.diff), as.factor(win)), fill=as.factor(win))) +
                 geom_boxplot()


      freq.plot <- ggplot(ml, aes(x=as.factor(variable), y=value,
                           group=interaction(variable, as.factor(winf)), fill=as.factor(winf))) +
                  geom_boxplot() +
                  facet_grid(A.geno~B.geno) +
                  ylab("Alt. allele freq") +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1))



      return(plot_grid(or.plot, freq.plot, ncol=1))

    }

    tmp <- plot.it()
    plot_grid(gwas.plot + theme(legend.position="none"), tmp, ncol=2)

  ### what genes?
    gff <- fread("/mnt/pricey_2/DaphniaGenomeAnnotation/Daphnia.aed.0.6.gff")

    tstrsplit(gff[V1==o[which.min(ppa)]$chr][V4>=o[which.min(ppa)]$start][V5<=o[which.min(ppa)]$stop][V3=="gene"]$V9, ";")[[1]]
    system("")


    library(SeqArray)
    genofile <- seqOpen("/mnt/pricey_2/DaphniaGenomeAnnotation/snpEff/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.annotated.gds")

    seqSetFilter(genofile, variant.id=m.inform[chr==o[which.min(ppa)]$chr][pos>=o[which.min(ppa)]$start][pos<=o[which.min(ppa)]$stop]$id)
    tmp <- seqGetData(genofile, "annotation/info/ANN")
    tmp.dt <- data.table(times=tmp$length, id=seqGetData(genofile, "variant.id"))
    tmp.l <- tmp.dt[,list(id=rep(id, times)), list(i=id)]
    tmp.l[,ann:=tmp$data]





















### genotype test

  ### test one
    m.l <- melt(m[use.chr==T,c("chr", "pos", "A.geno", "B.geno", "m.hat", "f.hat"),with=F],
                id.vars=c("chr", "pos", "A.geno", "B.geno"))


    ggplot(data=m.l[(A.geno==0 & B.geno==2) | (A.geno==2 & B.geno==0)],
            aes(x=as.factor(A.geno), y=value, group=interaction(variable, A.geno), color=variable)) +
    geom_boxplot() + ylab("Allele Frequency") + xlab("A's genotype")

  ### test two
    hist(m[A.geno==B.geno & m.hat!=0 & f.hat!=0]$mf_delta, breaks=100)
    hist(m[A.geno==0 & B.geno==2 & m.hat!=0 & f.hat!=0]$mf_delta, breaks=100)
    hist(m[A.geno==2 & B.geno==0 & m.hat!=0 & f.hat!=0]$mf_delta, breaks=100)
    hist(m[A.geno==0 & B.geno==1 & m.hat!=0 & f.hat!=0]$mf_delta, breaks=100)

### sliding window ancestry test
  chrs <- m[use.chr==T, list(min=min(pos), max=max(pos)), chr]
  window.bp <- 100000
  step.bp <- 500

  wins <- foreach(i=1:dim(chrs)[1], .combine="rbind")%dopar%{
    #i<-1
    print(i)
    tmp <- data.table(chr=chrs[i]$chr, start=seq(from=chrs[i]$min, to=chrs[i]$max-window.bp, by=step.bp))
    tmp[,stop:=start+window.bp]
    tmp
  }

  setkey(m, chr, pos)

  m.inform <- m[final.use==T][inform==T]

  o <- foreach(i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
    #i<-145167 ; i <- 58987
    print(paste(i, dim(wins)[1], sep=" / "))


    m.inform[,win:=0]
    m.inform[chr==wins$chr[i] & pos>=wins$start[i] & pos<=wins$stop[i], win:=1]

    if(plot.it==T) {
      diffPlot <- ggplot(data=m.inform, aes(x=as.factor(BA.geno.diff.polarize), y=mf_delta.polarize, group=interaction(BA.geno.diff.polarize, win), fill=as.factor(win))) + geom_boxplot()
      #ggplot(data=m.inform, aes(x=as.factor(BA.geno.diff), y=mf_delta, group=interaction(BA.geno.diff, win), fill=as.factor(win))) + geom_violin()

        # .polarize
      m.inform.l <- melt(m.inform[,c("chr", "pos", "win", "BA.geno.diff.polarize", "m.hat.polarize", "f.hat.polarize"), with=F],
      id.vars=c("chr", "pos", "win", "BA.geno.diff.polarize"))

      ggplot(data=m.inform.l, aes(x=as.factor(win), y=value, group=interaction(win, variable), fill=variable)) +
      geom_hline(yintercept=.5) + geom_boxplot() + facet_wrap(~BA.geno.diff.polarize)


      chrPlot.closeup <- ggplot(m.inform.l[win==1], aes(x=variable, y=as.numeric(as.factor(pos)), fill=as.factor(value))) + geom_tile() + coord_flip()

      chrPlot <- ggplot(m.inform.l[chr==unique(chr[win==1])], aes(x=variable, y=as.numeric(as.factor(pos)), fill=as.factor(value))) +
      geom_tile() + coord_flip() + geom_rect(aes(xmin=.5, xmax=2.5,
                                                 ymin=min(as.numeric(as.factor(pos))[win==1]),
                                                 ymax=max(as.numeric(as.factor(pos))[win==1])),colour="black", fill=NA, size=1)

       plot_grid(diffPlot, chrPlot.closeup, nrow=2)

    }

    mod <- summary(lm(mf_delta.polarize~BA.geno.diff.polarize*win, m.inform))


    data.table(win.i=i, chr=wins$chr[i], start=wins$start[i], stop=wins$stop[i], nSNPs=sum(m.inform$win),
               r2=mod$r.squared,
               beta.int=mod$coefficients[4,1],
               se=mod$coefficients[4,2],
               p=mod$coefficients[4,4])

    }

  o <- rbindlist(o)
  o[,pa:=p.adjust(p, "bonferroni")]
  save(o, file="~/fuckingRad.Rdata")

  load(file="~/fuckingRad.Rdata")


  ggplot(data=o[nSNPs>25], aes(y=t.diff, x=win.i)) + geom_point()

  ggplot(data=o[nSNPs>25], aes(y=-log10(t.p), x=win.i)) + geom_point()
  ggplot(data=o[nSNPs>25], aes(y=beta, x=win.i)) + geom_point()
  ggplot(data=o[nSNPs>25], aes(y=r2, x=win.i)) + geom_point()

  ggplot(data=o[nSNPs>25], aes(x=beta, y=-log10(p), color=nSNPs)) + geom_point()


  ggplot(data=o, aes(y=-log10(pa), x=win.i, color=chr)) + geom_point()


  #i<-14520 ; i <- 3854
  print(paste(i, dim(wins)[1], sep=" / "))

  tmp <-  m.l[J(data.table(chr=wins$chr[i], pos=wins$start[i]:wins$stop[i])), nomatch=0][A.geno!=B.geno]

  ggplot(data=tmp,
          aes(x=as.factor(A.geno), y=value, group=interaction(variable, A.geno), color=variable)) +
  geom_boxplot(position="dodge") + ylab("Allele Frequency") + xlab("A's genotype")






  plot_grid(A_AA.B_aa.plot , A_aa.B_AA.plot, labels=c("A:2; B:0", "A:0; B:2"))

  ggplot(data=m.l[A.geno==2 & B.geno==0], aes(x=interaction(A.geno, variable), y=value, color=variable)) + geom_boxplot()


  t1.male <- lm(value~A.geno+B.geno, m.l[variable=="m.hat"])
  t1.female <- lm(value~A.geno+B.geno, m.l[variable=="f.hat"])

  t0.mf <- lm(value~A.geno + B.geno + variable, m.l)
  t1.mf <- lm(value~A.geno*variable + B.geno*variable, m.l)



### likelihood
library(data.table)
library(foreach)
library(ggplot2)

  load("/mnt/sammas_storage/bergland-lab/alan/pool_AB.Rdata")
  m[!is.na(A),A.geno := unlist(sapply(m[!is.na(A)]$A, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[!is.na(B),B.geno := unlist(sapply(m[!is.na(B)]$B, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[,BA.geno.diff := B.geno - A.geno]
  m[,inform:=!((A.geno==0 & B.geno==0) | (A.geno==2 & B.geno==2))]

  m[!is.na(A.geno) & !is.na(B.geno), AB.genoclass := paste(A.geno, B.geno, sep="_")]

  exp.tab <- data.table(exp.freq=c(0, .25, .5, .25, .5, .75, .5, .75, 1),
                        AB.genoclass=c("0_0", "0_1", "0_2", "1_0", "1_1", "1_2", "2_0", "2_1", "2_2"),
                        sign=c(0, 1, 1, 0, 0, 1, 0, 0, 0))
  setkey(m, AB.genoclass)
  setkey(exp.tab, AB.genoclass)

  m <- merge(m, exp.tab)
  m <- m[exp.freq!=0 & exp.freq!=1]

  ggplot() +
  geom_boxplot(data=m[effRD_D8Male1>40][final.use==T],
               aes(x=as.factor(AB.genoclass), y=effPA_D8Male1)) +
  geom_point(data=exp.tab,
             aes(as.factor(AB.genoclass), y=exp.freq), color="red")


  fun <- function(c.i, samp="D8Male1", tmp) {

      ll <- sum(dbinom(as.matrix(tmp[,paste("effRD_", samp, sep=""), with=F] * tmp[,paste("effPA_", samp, sep=""), with=F])[,1],
                     as.matrix(tmp[,paste("effRD_", samp, sep=""), with=F])[,1],
                     plogis(qlogis(tmp$exp.freq)+(tmp$sign - c.i)), log=T),
                 na.rm=T)

    #data.table(ll=ll, c=c.i)
    -ll
  }


  o <- foreach(c.i=seq(from=-1, to=1, by=.05), .combine="rbind")%dopar%{
    fun(samp="D8Male2", c=c.i)
  }
  plot(ll~c, o)
  o[which.max(ll)]


  o <- nlm(fun, p=.4, hessian=TRUE)


  chrs <- m[use.chr==T, list(min=min(pos), max=max(pos)), chr]
  window.bp <- 100000
  step.bp <- 10000

  wins <- foreach(i=1:dim(chrs)[1], .combine="rbind")%dopar%{
    #i<-1
    print(i)
    tmp <- data.table(chr=chrs[i]$chr, start=seq(from=chrs[i]$min, to=chrs[i]$max-window.bp, by=step.bp))
    tmp[,stop:=start+window.bp]
    tmp
  }
  wins[,i:=1:dim(wins)[1]]

  setkey(m, chr, pos)
  o <- foreach(i=c(1:dim(wins)[1]), .errorhandling="remove")%dopar%{
    print(i)
    tmp.dat <- m[J(data.table(chr=wins$chr[i], pos=wins$start[i]:wins$stop[i])), nomatch=0]

    foreach(x=c("D8Male1", "D8Male2", "D8PE1", "D8PE2"), .combine="rbind")%do%{
      o <- nlm(fun, p=.4, tmp=tmp.dat, samp=x)
      data.table(i=i, c=o$estimate, n=dim(tmp.dat)[1], samp=x)
    }

  }
  o <- rbindlist(o)
  o <- merge(o, wins, by="i")


  ggplot(data=o[n>50], aes(x=i, y=c, group=samp, color=as.factor(grepl("Male", samp)))) + geom_line()

  o.ag <- o[n>50,list(diff=(c[samp=="D8PE2"] -
                            c[samp=="D8PE1"])), list(i, chr, start, stop)]

  plot(diff~i, o.ag)

#############
### PLOTS ###
#############

  ### AB & poolMF comparison.
    ggplot(data=m.ag) +
    geom_line(aes(x=delta, y=log2(or))) +
    geom_line(aes(x=delta, y=log2(uci)), linetype="dashed") +
    geom_line(aes(x=delta, y=log2(lci)), linetype="dashed")

  ### focusing in on regions
    plotRegion <- function(chr, start, stop) {
      #chr <- "Scaffold_7757_HRSCAF_8726"; start <- 8244208; stop <- 8576809

      tmp <- m[chr==chr][pos>=start][pos<=stop]

      ggplot(data=tmp[BA_delta!=0], aes(x=pos, y=abs(mf_delta))) + geom_line()
      ggplot(data=tmp[BA_delta!=0], aes(x=(BA_delta), y=(mf_delta))) + geom_point()

    }
