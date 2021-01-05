### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(ggtern)
  library(viridis)

### load output data from HWE_simulations.scatter.rivanna.R
  fl <- system("ls /scratch/aob2x/daphnia_hwe_sims/hweSimOut/*.csv", inter=T)

  hwe.ag.m <- foreach(i=fl)%do%{
    print(i)
    fread(i)
  }
  hwe.ag.m <- rbindlist(hwe.ag.m)

  hwe.ag.m.ag <- hwe.ag.m[,list(diff.mu=mean(diff), diff.sd=sd(diff), n.obs=mean(n.obs), n.sim=mean(n.sim)),
                           list(fAA, fAa, faa, py)]

  save(hwe.ag.m.ag, hwe.ag.m, file="/nv/vol186/bergland-lab/Daphnia_HWE/hwe.ag.m.ag")













### OLD & DEFUNCT?  




  sum(hwe.ag.m.ag[py=="D8.2016.2017.2018.2019"]$diff.mu) / sum(hwe.ag.m.ag[py=="D8.2016.2017.2018.2019"]$n.obs)

  hist(hwe.ag.m.ag[py=="D8.2016.2017.2018.2019"]$diff.mu)



tmp <- hwe.ag.m[py=="D8.2016.2017.2018.2019"]
tmp[which.max(diff)]

ggplot() +
geom_point(data=hwe.ag.m[py=="D8.2016.2017.2018.2019"][order(diff, decreasing=F)][diff!=0],
            aes(x=fAA, y=fAa, z=faa, color=n.obs)) +
coord_tern(expand=T) + limit_tern(T = 1.05, L = 1.05, R = 1.05) +
scale_color_viridis() + scale_fill_viridis()




ggplot() +
geom_point(data=hwe.ag.m[py=="D8.DBunk.2016.2017.2018.2019"][order(diff, decreasing=F)][diff!=0],
            aes(x=fAA, y=fAa, z=faa, color=sign(diff)*(abs(diff)))) +
coord_tern(expand=T) + limit_tern(T = 1.05, L = 1.05, R = 1.05) +
scale_color_viridis() + scale_fill_viridis()



tmp <- hwe.ag.m[py=="D8.2016.2017.2018.2019"]




table(hwe.ag$pond, hwe.ag$year)
hwe.ag[,py:=paste(pond, year, sep=".")]

py.1 <- ggplot() +
        geom_point(data=hwe.ag[py=="D8.2016.2017.2018.2019"][order(n, decreasing=F)], aes(x=fAA, y=fAa, z=faa, color=(n))) +
        coord_tern(expand=T) + limit_tern(T = 1.05, L = 1.05, R = 1.05) +
        scale_color_viridis() + scale_fill_viridis()

py.2 <- ggplot() +
        geom_point(data=hwe.ag[py=="DBunk.2017"][order(n, decreasing=F)], aes(x=fAA, y=fAa, z=faa, color=(n))) +
        coord_tern(expand=T) + limit_tern(T = 1.05, L = 1.05, R = 1.05) +
        scale_color_viridis() + scale_fill_viridis()

        ggplot() +
        geom_point(data=hwe.ag[py=="D8.DBunk.2016.2017.2018.2019"][order(n, decreasing=F)], aes(x=fAA, y=fAa, z=faa, color=(n))) +
        coord_tern(expand=T) + limit_tern(T = 1.05, L = 1.05, R = 1.05) +
        scale_color_viridis() + scale_fill_viridis()


plot_grid(py.1, py.2, labels=c("D8 / 2017.2018", "DBunk / 2017"))


### is there heterogeneity among chromosomes? YES.
hwe.stat[,py:=paste(pond, year, sep=".")]

hwe.stat[,odd:=NA]

for(py.i in c("DBunk.2017", "D8.2017.2018")) {
  print(py.i)
  ideal <- hwe.ag[py==py.i][fAA!=0.5][which.max(n)]$fAA
  buffer <- 0.025

  hwe.stat[py==py.i, odd:=F]
  hwe.stat[py==py.i & fAA>=(ideal-buffer) & fAA<=(ideal+buffer) & fAa>=(1-ideal-buffer) & fAa<=(1-ideal+buffer) & faa<=buffer, odd:=T]
}

chisq.test(table(hwe.stat[py=="DBunk.2017"]$chr, hwe.stat[py=="DBunk.2017"]$odd))
chisq.test(table(hwe.stat[py=="D8.2017.2018"]$chr, hwe.stat[py=="D8.2017.2018"]$odd))


hwe.chr <- foreach(py.i=c("DBunk.2017", "D8.2017.2018"), .combine="rbind")%do%{
  temp <- hwe.stat[py==py.i,list(n.ZW=sum(odd==T), n.auto=sum(odd==F)), list(chr, py)]
  temp[,exp.n.ZW:=mean(hwe.stat[py==py.i]$odd==T) * (n.ZW+n.auto)]
  temp[,en:=(n.ZW - exp.n.ZW)/exp.n.ZW]
  temp
}

plot(hwe.chr[py=="DBunk.2017"]$en ~ hwe.chr[py=="D8.2017.2018"]$en)

### sliding window
hwe.stat <- na.omit(hwe.stat)

window.bp <- 100000
step.bp <- window.bp/20
setkey(hwe.stat, chr)

wins <- foreach(chr.i=unique(hwe.stat$chr), .combine="rbind")%do%{
  data.table(start=seq(from=min(hwe.stat[J(chr.i)]$pos),
                               to=max(hwe.stat[J(chr.i)]$pos) - window.bp,
                               by=step.bp),
                     stop=seq(from=min(hwe.stat[J(chr.i)]$pos),
                              to=max(hwe.stat[J(chr.i)]$pos) - window.bp,
                              by=step.bp)+window.bp,
                      chr=chr.i)
}

setkey(hwe.stat, chr, pos)

o <- foreach(i=1:dim(wins)[1])%dopar%{

  print(paste(i, dim(wins)[1], sep=" / "))

  temp <- hwe.stat[J(data.table(chr=wins[i]$chr, pos=c(wins[i]$start:wins[i]$stop))), nomatch=0]

  temp[,list(rate=mean(odd==T), n=length(odd),
             i=i, ch=wins[i]$chr, min.pos=wins[i]$start, max.pos=wins[i]$stop),
         py]
}

o <- rbindlist(o)

o.br <- hwe.stat[,list(base.rate=mean(odd==T)), list(py)]
setkey(o, py)
setkey(o.br, py)
o <- merge(o, o.br)
o[,en := log2(rate/base.rate)]
#  o[,en := (rate-base.rate)/base.rate]

ow <- dcast(o[n>50], i~py, value.var="en")
plot(I(2^(D8.2017.2018))~I(2^(DBunk.2017)), ow)
abline(0,1)

ow.fet <- foreach(i=seq(from=1, to=10, by=.5), .combine="rbind", .errorhandling="remove")%dopar%{
  foreach(j=seq(from=1, to=10, by=.5), .combine="rbind", .errorhandling="remove")%do%{

    print(paste(i, i, sep=" . "))
    fet <- fisher.test(2^ow$DBunk.2017>i, 2^ow$D8.2017.2018>j)
    data.table(i=i, j=j, or=fet$estimate, p=fet$p.value)
  }
}
ggplot(data=ow.fet, aes(x=i, y=j, fill=(or))) + geom_tile() + scale_color_viridis() + scale_fill_viridis()


ggplot(data=o[n>50], aes(x=i, y=rate, color=ch)) + facet_grid(py~.) +
geom_vline(xintercept=ow[2^D8.2017.2018>10 & 2^DBunk.2017>10]$i) +
geom_point()







### what is the peak region?
  ggplot(data=hwe[J(data.table(chr=o[n>50][which.max(en)]$ch ,
                               pos=o[n>50][which.max(en)]$min.pos:o[n>50][which.max(en)]$max.pos)), nomatch=0],
          aes(x=fAA, y=fAa, z=faa)) + coord_tern() + geom_point()




hwe <- hwe[p<.99]

ggplot() +
geom_point(data=hwe[p>1e-4][order(f)], aes(x=fAA, y=fAa, z=faa, color=f), size=1.5) +
geom_point(data=hwe[p<1e-3][f>0], aes(x=fAA, y=fAa, z=faa), color="yellow", size=.05) +
geom_point(data=hwe[p<1e-3][f<0], aes(x=fAA, y=fAa, z=faa), color="red", size=.05) +
coord_tern() +
limit_tern(T = 0.7, L = 0.7, R = 0.7)






Tlim=c(-10, 110), Llim=c(-.1, 1.1)*100, Rlim=c(-.1, 1.1)*100

setkey(hwe, id)
setkey(snp.dt, id)
hwe <- merge(hwe, snp.dt, by.x="variant.id", by.y="id")



### sliding window
setkey(clones.af, chr)


window.bp <- 100000
step.bp <- 10000

o <- foreach(chr.i=unique(clones.af$chr), .combine="rbind", .errorhandling="remove")%do%{
  ### test
    #chr.i <- "Scaffold_1863_HRSCAF_2081"
    chr.i <- "Scaffold_9201_HRSCAF_10758"

  ### define windows
    wins <- data.table(start=seq(from=min(clones.af[J(chr.i)]$pos),
                                 to=max(clones.af[J(chr.i)]$pos) - window.bp,
                                 by=step.bp),
                       stop=seq(from=min(clones.af[J(chr.i)]$pos),
                                to=max(clones.af[J(chr.i)]$pos) - window.bp,
                                by=step.bp)+window.bp)

  ### iterate through windows
    o <- foreach(k=1:dim(wins)[1], .combine="rbind", .errorhandling="remove")%dopar%{
      print(paste(chr.i, k, dim(wins)[1], sep=" / "))


      ### extract out SNP IDs in window
        snp.temp <- clones.af[J(chr.i)][pos>=wins[k]$start & pos<wins[k]$stop]$id

      ### LD stuff
      #temp <- snpgdsLDMat(genofile,
      #            sample.id=clones[year==2017]$clone,
      #            snp.id=snp.temp,
      #            method="r", verbose=F, num.thread=1)

      #LD.expand <- as.data.table(expand.grid(temp$LD, KEEP.OUT.ATTRS = FALSE))
      #LD.expand[,id1:=rep(snp.temp, length(snp.temp))]
      #LD.expand[,id2:=rep(snp.temp, each=length(snp.temp))]

      #LD.expand <- merge(LD.expand, clones.af[,c("id", "pos"), with=F], by.x="id1", by.y="id")
      #LD.expand <- merge(LD.expand, clones.af[,c("id", "pos"), with=F], by.x="id2", by.y="id")
      #LD.expand[,dist := abs(pos.x - pos.y)]

      ### HWE stuff
        #hwe.p <- snpgdsHWE(genofile,
        #            sample.id=clones[year==2017]$clone,
        #            snp.id=snp.temp)

        seqSetFilter(genofile, sample.id=clones[year=="2017"]$clone, variant.id=snp.temp)

        #f <- inbreedCoeff(genofile, margin="by.variant")
        hwe <- as.data.table(hwe(genofile))

      ### format output

      data.table(start=wins[k]$start, stop=wins[k]$stop, chr=chr.i,
                 #LD.median=median(temp$LD^2, na.rm=T),
                 #LD.mean=mean(temp$LD^2, na.rm=T),
                 #LD.decay.beta=summary(lm(Var1~dist, LD.expand))$coef[2,1],
                 #LD.decay.p=summary(lm(Var1~dist, LD.expand))$coef[2,4],
                 HWE.med.p = median(hwe[p<.999]$p),
                 f = median(hwe[p<.999]$f),

                 nSNPs=length(snp.temp),
                 win.id=k)


    }
    o[,med:=start/2 + stop/2]

  ### return
    return(o)
}




ggplot(data=o[nSNPs>15], aes(x=med, y=-log10(HWE.med.p), color=chr)) + geom_point() + facet_wrap(~chr, scale="free_x")
ggplot(data=o[nSNPs>15], aes(x=med, y=f, color=chr)) + geom_point() + facet_wrap(~chr, scale="free_x")





#### SW averages


window.bp <- 50000
step.bp <- 5000

o <- foreach(chr.i=unique(clones.af$chr), .combine="rbind", .errorhandling="remove")%do%{
### define windows
  wins <- data.table(start=seq(from=min(clones.af[J(chr.i)]$pos),
                               to=max(clones.af[J(chr.i)]$pos) - window.bp,
                               by=step.bp),
                     stop=seq(from=min(clones.af[J(chr.i)]$pos),
                              to=max(clones.af[J(chr.i)]$pos) - window.bp,
                              by=step.bp)+window.bp)

### iterate through windows
  o <- foreach(k=1:dim(wins)[1], .combine="rbind", .errorhandling="remove")%dopar%{
    print(paste(chr.i, k, dim(wins)[1], sep=" / "))

    snp.temp <- clones.af[J(chr.i)][pos>=wins[k]$start & pos<wins[k]$stop]$id
    hwe.temp <- hwe[J(snp.temp)]

    hwe.temp <- hwe.temp[nAA>10 | naa>10]

    ### format output

    data.table(start=wins[k]$start, stop=wins[k]$stop, chr=chr.i,
               start.id=min(snp.temp), stop.id=max(snp.temp),
               #LD.median=median(temp$LD^2, na.rm=T),
               #LD.mean=mean(temp$LD^2, na.rm=T),
               #LD.decay.beta=summary(lm(Var1~dist, LD.expand))$coef[2,1],
               #LD.decay.p=summary(lm(Var1~dist, LD.expand))$coef[2,4],
               HWE.med.p = median(hwe.temp[p<.999]$p),
               f = median(hwe.temp[p<.999]$f),
               nAA.mu=mean(hwe.temp$nAA),
               nAa.mu=mean(hwe.temp$nAa),
               naa.mu=mean(hwe.temp$naa),
               nSNPs=length(snp.temp),
               win.id=k)


  }
  o[,med:=start/2 + stop/2]

### return
  return(o)
}

ggplot() +
geom_point(data=o, aes(x=med, y=f, color=chr)) +
geom_point(data=o[HWE.med.p<1e-3], aes(x=med, y=f), size=.1) +
facet_wrap(~chr, scale="free_x") + geom_hline(yintercept=0)

ggplot(data=o, aes(x=med, y=-log10(HWE.med.p), color=chr)) + geom_point() + facet_wrap(~chr, scale="free_x")












hwe[chr=="Scaffold_6786_HRSCAF_7541"][pos>=5933032 & pos<=6283032]


  plot(f~afreq, hwe[(nAA > 5 | naa>5) & nAa>5])
