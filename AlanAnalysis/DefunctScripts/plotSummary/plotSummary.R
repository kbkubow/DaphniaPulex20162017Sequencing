### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R


  ### libraries
    library(data.table)
    library(foreach)
    library(cowplot)
    library(SeqArray)

  ### load output from harp (generated by `plotHarp.R`)
    #load(file="/scratch/aob2x/daphnia_hwe_sims/harp_pools/summarizedOut/harpWins.Rdata")
    load("/mnt/sammas_storage/bergland-lab/alan/harpWins.Rdata")
    setkey(o, chr)
    harpWins <- o
    rm(o)

  ### load sliding window plot (from above)
    load(file="~/slidingWindow_oddsRatio.Rdata")

  ### load precomputed mInform file (from above)
    load(file="~/mInform.Rdata")

  ### open connection to annotated GDS file
    #genofile <- seqOpen("/mnt/pricey_2/DaphniaGenomeAnnotation/snpEff/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.annotated.gds")

  ### load haplotype priors
    dat <- fread("/mnt/sammas_storage/bergland-lab/alan/testTrio.consensus.header.phase.csv")

    setnames(dat, c(1,2), c("chr", "pos"))

    dat <- dat[!(A=="1/1" & B=="1/1")][!(A=="0/0" & B=="0/0")]

    dat.p <- dat[,list(ref=REF, alt=ALT,
                        A1=as.numeric(substr(A, 0, 1)),
                        A2=as.numeric(substr(A, 3, 3)),
                        B1=as.numeric(substr(B, 0, 1)),
                        B2=as.numeric(substr(B, 3, 3))),
                    list(chr, pos)]

  ### load pulicaria
    #puli <- fread("/mnt/sammas_storage/bergland-lab/alan/pulicaria.sort.D84a.rmDup.aseReadCounter.allvariant.delim")
    #setnames(puli, c("contig", "position"), c("chr", "pos"))
    #setkey(puli, chr, pos)
    #puli[,puli.geno:=round(refCount/totalCount*10)/10]
    #puli[puli.geno!=1 & puli.geno!=0, puli.geno:=.5]

    #### playing
      #setkey(dat.p, chr, pos)
      #foo <- merge(puli[,c("chr", "pos", "puli.geno")], dat.p, all.x=T, all.y=T)
      #prop.table(table(foo$puli.geno==foo$A.geno))
      #prop.table(table(foo$puli.geno==foo$B.geno))

      #foo.l <- melt(foo, id.vars=c("chr", "pos", "puli.geno", "ref", "alt"))
      #foo.l.ag <- foo.l[puli.geno%in%c(0,1), list(IBS=mean(value==puli.geno, na.rm=T), .N), list(chr, variable)]

      #ggplot(data=foo.l.ag[N>1000], aes(x=variable, y=IBS, color=chr)) + geom_point() + coord_flip()

    load("/mnt/sammas_storage/bergland-lab/alan/puAg.Rdata")
    setkey(pu.ag, chr, pos)

    #### playing
      setkey(dat.p, chr, pos)
      foo <- merge(pu.ag[,c("chr", "pos", "puli.geno")], dat.p)
      foo[,puli.geno:=puli.geno/2]
      prop.table(table(foo$puli.geno==foo$A.geno))
      prop.table(table(foo$puli.geno==foo$B.geno))

      foo.l <- melt(foo, id.vars=c("chr", "pos", "puli.geno", "ref", "alt"))
      foo.l.ag <- foo.l[puli.geno%in%c(0,1), list(IBS=mean(value==puli.geno, na.rm=T), .N), list(chr, variable)]

      ggplot(data=foo.l.ag[N>1000], aes(x=variable, y=IBS, color=chr)) + geom_point() + coord_flip()



  ### generte windows
    setkey(dat.p, chr, pos)
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


  ### polarized Manhattan plot
    o[,ppa:=(sign(beta.int)*log10(pa))]
    gwas.plot <- ggplot(data=o[nSNPs>25], aes(y=(ppa), x=win.i, color=chr)) + geom_line()

  ### plotting functions
    plot.or <- function(i=o[which.min(ppa)]$win.i) {
      #i <- wins[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(start-8660157))]$i
      #i <- o[nSNPs>25][which.max(ppa)]$win.i
      #i <- o[win.i>5000][which.min(ppa)]$win.i

      m.inform[,win:=0]
      m.inform[chr==wins$chr[i] & pos>=wins$start[i] & pos<=wins$stop[i], win:=1]

      ml <- melt(m.inform[,c("chr", "pos", "win", "D8Male.f", "D8PE.f", "BA.geno.diff", "A.geno", "B.geno"), with=F],
                 id.vars=c("chr", "pos", "win", "BA.geno.diff", "A.geno", "B.geno"))
      ml[,winf:=c("genome", "targetWindow")[win+1]]
      ml[,exp_hom:= 1 - 2*value*(1-value)]

      ggplot(m.inform, aes(x=as.factor((BA.geno.diff)), y=log2(or.fold),
                           group=interaction((BA.geno.diff), as.factor(win)), fill=as.factor(win))) +
                 geom_boxplot()
               }

    plot.freq <- function(i=o[which.min(ppa)]$win.i) {
     #i <- wins[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(start-8660157))]$i
     #i <- o[nSNPs>25][which.max(ppa)]$win.i
     #i <- o[win.i>5000][which.min(ppa)]$win.i

     m.inform[,win:=0]
     #m.inform[chr==wins$chr[i] & pos>=wins$start[i] & pos<=wins$stop[i], win:=1]
     m.inform[chr==o[win.i==i]$chr & pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop, win:=1]



     ml <- melt(m.inform[,c("chr", "pos", "win", "D8Male.f", "D8PE.f", "BA.geno.diff", "A.geno", "B.geno"), with=F],
                id.vars=c("chr", "pos", "win", "BA.geno.diff", "A.geno", "B.geno"))
     ml[,winf:=c("genome", "targetWindow")[win+1]]
     ml[,exp_hom:= 1 - 2*value*(1-value)]


     ggplot(ml, aes(x=as.factor(variable), y=value,
                           group=interaction(variable, as.factor(winf)), fill=as.factor(winf))) +
                  geom_boxplot() +
                  facet_grid(A.geno~B.geno) +
                  ylab("Alt. allele freq") +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1))
                }

    plot.hap <- function(i=o[which.min(ppa)]$win.i) {

     ### haplotype painting plot
      ### extract out A&B haplotypes
       tmp <- dat.p[chr==o[win.i==i]$chr & pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]

       tmp <- merge(tmp, pu.ag[,c("chr", "pos", "puli.geno"), with=F], all.x=T)
       tmp[,puli.geno:=puli.geno/2]

       tmp.ag <- tmp[,list(i=seq(from=min(pos), to=max(pos), length.out=length(pos)), pos=pos), list(chr)]
       setkey(tmp.ag, chr, pos)
       setkey(tmp, chr, pos)
       tmp <- merge(tmp, tmp.ag)



       tmp[A1==0, A1.allele:=alt]
       tmp[A1==1, A1.allele:=ref]
       tmp[A2==0, A2.allele:=alt]
       tmp[A2==1, A2.allele:=ref]

       tmp[B1==0, B1.allele:=alt]
       tmp[B1==1, B1.allele:=ref]
       tmp[B2==0, B2.allele:=alt]
       tmp[B2==1, B2.allele:=ref]

       tmp[puli.geno==0, puli.allele:=alt]
       tmp[puli.geno==1, puli.allele:=ref]
       tmp[puli.geno==.5, puli.allele:=ref]

       tmp2 <- tmp
       tmp2[tmp$B2==0, A1:=abs(1-A1)]
       tmp2[tmp$B2==0, A2:=abs(1-A2)]
       tmp2[tmp$B2==0, B1:=abs(1-B1)]
       tmp2[tmp$B2==0, B2:=abs(1-B2)]
       tmp2[tmp$B2==0, puli.geno:=abs(1-puli.geno)]



       tmp2[tmp$B2==0 & A1==0, A1.allele:=ref]
       tmp2[tmp$B2==0 & A1==1, A1.allele:=alt]
       tmp2[tmp$B2==0 & A2==0, A2.allele:=ref]
       tmp2[tmp$B2==0 & A2==1, A2.allele:=alt]

       tmp2[tmp$B2==0 & B1==0, B1.allele:=ref]
       tmp2[tmp$B2==0 & B1==1, B1.allele:=alt]
       tmp2[tmp$B2==0 & B2==0, B2.allele:=ref]
       tmp2[tmp$B2==0 & B2==1, B2.allele:=alt]

       tmp2[tmp$B2==0 & puli.geno==0, puli.allele:=ref]
       tmp2[tmp$B2==0 & puli.geno==1, puli.allele:=alt]
       tmp2[tmp$B2==0 & puli.geno==0.5, puli.allele:=ref]



       tmpl <- melt(tmp2, id.vars=c("chr", "pos", "ref", "alt", "i"),
                      measure=list(c("A1", "A2", "B1", "B2", "puli.geno"),
                                   c("A1.allele", "A2.allele", "B1.allele", "B2.allele", "puli.allele")),
                       value.name=c("value", "allele"))
       tmpl[,variable:=c("A1", "A2", "B1", "B2", "puli.geno")[variable]]

     ### fold in annotations
       seqSetFilter(genofile, variant.id=m.inform[chr==o[win.i==i]$chr][pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]$id)
       ann.tmp <- seqGetData(genofile, "annotation/info/ANN")
       ann.tmp.dt <- data.table(times=ann.tmp$length, id=seqGetData(genofile, "variant.id"))
       ann.tmp.l <- ann.tmp.dt[,list(id=rep(id, times)), list(i=id)]
       ann.tmp.l[,ann:=ann.tmp$data]

       classes <- c("stop|mis")
       ann.tmp.l <- ann.tmp.l[grepl(classes, ann)]

       setkey(ann.tmp.l, id)
       setkey(m.inform, id)

       ann.tmp.l[,allele:=tstrsplit(ann, "\\|")[[1]]]
       ann.tmp.l[,mut.class:=tstrsplit(ann, "\\|")[[2]]]
       ann.tmp.l[,gene:=tstrsplit(ann, "\\|")[[4]]]
       ann.tmp.l <- ann.tmp.l[,list(.N), list(id, allele, mut.class, gene)]
       ann.tmp.m <- merge(ann.tmp.l[,c("id","allele", "mut.class")], m.inform[,c("chr", "pos", "id", "A.geno", "B.geno")], by="id")

       setkey(tmpl, chr, pos, allele)
       setkey(ann.tmp.m, chr, pos, allele)

       foo <- merge(tmpl, ann.tmp.m, all.x=T)

       haps <- ggplot() +
                geom_tile(data=foo, aes(x=variable, y=i, fill=as.factor(abs(1-value)))) +
                geom_jitter(data=foo[!is.na(mut.class)], aes(x=variable, y=i, shape=mut.class, size=mut.class), color="black", width=.1) +
                coord_flip()  +
                theme(legend.position = "none")

     ### annotations
      seqSetFilter(genofile, variant.id=m.inform[chr==o[win.i==i]$chr][pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]$id)
      ann.tmp <- seqGetData(genofile, "annotation/info/ANN")
      ann.tmp.dt <- data.table(times=ann.tmp$length, id=seqGetData(genofile, "variant.id"))
      ann.tmp.l <- ann.tmp.dt[,list(id=rep(id, times)), list(i=id)]
      ann.tmp.l[,ann:=ann.tmp$data]

      setkey(ann.tmp.l, id)
      setkey(m.inform, id)

      ann.tmp.m <- merge(ann.tmp.l, m.inform[,c("chr", "pos", "id", "A.geno", "B.geno")])

      ann.tmp.m[grepl("stop", ann)]$pos[1]
      ann.tmp.m[grepl("missense", ann)]

     ### haplotype frequency plot
      pad <- 0
      o[,mid:=start/2 + stop/2]
      hapFreq.tmp <- harpWins

      hapFreq.tmp[,targetWin:=0]
      hapFreq.tmp <- hapFreq.tmp[chr==o[win.i==i]$chr & start<=(o[win.i==i]$mid-pad) & stop>=(o[win.i==i]$mid+pad), targetWin:=1]

      hapFreq.tmp.l <- melt(hapFreq.tmp, id.vars=c("chr", "start", "stop", "pool", "targetWin"))
      hapFreq.tmp.l[grepl("Male", pool),pool.f:="D8Male.f"]
      hapFreq.tmp.l[grepl("PE", pool),pool.f:="D8PE.f"]

      hapFreq.tmp.lag.ag <- hapFreq.tmp.l[,list(freq=mean(value, na.rm=T), sd=sd(value, na.rm=T)), list(allele=variable, pool=factor(pool, levels=rev(c("D8PE1", "D8PE2", "D8Male1", "D8Male2"))), targetWin)]

      haplofreq.plot <- ggplot() +
      geom_bar(data=hapFreq.tmp.lag.ag[targetWin==1], aes(x=allele, y=freq, group=pool, fill=pool), stat="identity", position="dodge") +
      geom_point(data=hapFreq.tmp.lag.ag[targetWin==0], aes(x=allele, y=freq, group=pool), position=position_dodge(1)) +
      geom_errorbar(data=hapFreq.tmp.lag.ag[targetWin==0], aes(x=allele, ymin=freq-sd, ymax=freq+sd, group=pool), position=position_dodge(1)) +
      coord_flip() + guides(fill = guide_legend(reverse = TRUE))

      ### map
        map <- ggplot() +
        geom_point(data=tmp, aes(y=1, x=seq(from=min(pos), to=max(pos), length.out=length(pos)))) +
        geom_point(data=tmp, aes(y=0, x=pos)) +
        geom_segment(data=tmp, aes(y=0, yend=1, x=pos, xend=seq(from=min(pos), to=max(pos), length.out=length(pos))), size=.01) +
        scale_size_manual()

      ### density
        dens <- ggplot() +
                geom_histogram(data=o, aes(nSNPs)) +
                geom_vline(data=o[win.i==i], aes(xintercept=nSNPs), color="red")


     plot_grid(haps, haplofreq.plot, map, dens, ncol=2)

    }

    plot.hap.freq <- function(i=o[which.min(ppa)]$win.i, bi, aj) {
      message(paste(i, bi, aj, sep=" / "))
      #i <- o[which.min(ppa)]$win.i
      #
      #x11()
      #bi<-1; aj<-0


      #rm(pooled.tmp, haps.tmp, pool.hap)

      pooled.exp <- m.inform[,list(freq.mu=c(mean(D8Male.f, na.rm=T), mean(D8PE.f, na.rm=T)), variable=c("D8Male.f", "D8PE.f")), list(A.geno, B.geno)]

      pooled.tmp <- m.inform[chr==o[win.i==i]$chr & pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]
      haps.tmp <- dat.p[chr==o[win.i==i]$chr & pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]


      setkey(pooled.tmp, chr, pos)
      setkey(haps.tmp, chr, pos)
      pool.hap <- merge(pooled.tmp[,c("chr", "pos", "A.geno", "B.geno", "D8Male.f", "D8PE.f"), with=F], haps.tmp, all=T)

      pool.hap[,A.haploSum:=abs(2-(A1+A2))]
      pool.hap[,B.haploSum:=abs(2-(B1+B2))]

      print(table(A=pool.hap$A.geno, B=pool.hap$B.geno))

      ### check
        #table(pool.hap$A.geno==pool.hap$A.haploSum); table(pool.hap$B.geno==pool.hap$B.haploSum)

      Bi.Aj <- na.omit(pool.hap[B.geno==bi][A.geno==aj])
      message(dim(Bi.Aj)[1])

      if(dim(Bi.Aj)[1]!=0) {
        Bi.Aj.freq.l <- melt(Bi.Aj, id.vars=c("chr", "pos"), measure.vars=c("D8Male.f", "D8PE.f"))
        Bi.Aj.haps.l <- melt(Bi.Aj, id.vars=c("chr", "pos"), measure.vars=c("A1", "A2", "B1", "B2"))

        freq.plot <- ggplot() +
                     geom_boxplot(data=Bi.Aj.freq.l, aes(x=as.factor(variable), y=value,
                                           group=variable)) +
                     ylab("Alt. allele freq") +
                     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                     ggtitle(paste("A:", aj, " B:", bi, "   nSNPs: ", dim(Bi.Aj)[1], sep="")) +
                     ylim(0,1) +
                     geom_point(data=pooled.exp[A.geno==aj & B.geno==bi], aes(x=as.factor(variable), y=freq.mu), color="red")

        hap.plot <- ggplot(data=Bi.Aj.haps.l, aes(x=variable, y=i, fill=as.factor(abs(1-value)))) + geom_tile() + coord_flip()  + theme(legend.position = "top")
        return(plot_grid(freq.plot, hap.plot))
      } else {
        freq.plot <- ggplot(data=data.table(x=1, y=1), aes(x=x, y=x))
        hap.plot <- freq.plot
        return(plot_grid(freq.plot, hap.plot))
      }

    }

    ii <- o[which.max(ppa)]$win.i
    ii <- o[win.i<5000][which.min(ppa)]$win.i
    ii <- o[which.min(ppa)]$win.i
  #  ii <- sample(o$win.i, 1)

    b1a0 <- plot.hap.freq(i=ii, bi=1, aj=0)
    b0a1 <- plot.hap.freq(i=ii, bi=0, aj=1)
    b2a1 <- plot.hap.freq(i=ii, bi=2, aj=1)
    b1a2 <- plot.hap.freq(i=ii, bi=1, aj=2)
    b1a1 <- plot.hap.freq(i=ii, bi=1, aj=1)

    x11(); plot_grid(b1a0, b0a1, b2a1, b1a2, ncol=2)

    hist(o$nSNPs)
    abline(v=o[which.min(ppa)]$nSNPs, col="red")



    ii <- o[which.max(ppa)]$win.i
    ii <- o[win.i<5000][which.min(ppa)]$win.i
    ii <- o[which.min(ppa)]$win.i

    plot.hap(i=ii)



    hapPlot <- plot.hap(i=o[which.min(ppa)]$win.i)
    freqPlot <- plot.freq(i=o[which.min(ppa)]$win.i)

    plot_grid(freqPlot, hapPlot, ncol=2, rel_widths=c(2,1))

   plot.hap(i=o[win.i<5000][which.min(ppa)]$win.i)
   plot.hap(i=o[which.max(ppa)]$win.i)

   plot.hap(i=sample(o$win.i, 1))
   plot.hap(i=o[which.max(nSNPs)]$win.i)
