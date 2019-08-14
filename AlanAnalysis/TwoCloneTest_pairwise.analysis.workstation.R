### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(ggtern)
  library(viridis)
  library(SeqArray)
  library(regioneR)
  library(doMC)
  registerDoMC(20)

### load HWE summary data (HWE_simulations.gather.rivanna.R makes this file)
  load(file="/mnt/sammas_storage/bergland-lab/Daphnia_HWE/hwe.ag.m.ag")

### open HWE & genotype dist statistic data
  ### cp /scratch/aob2x/daphnia_hwe_sims/hwe_stat.Rdata /nv/vol186/bergland-lab/Daphnia_HWE/hwe_stat.Rdata
  load(file="/mnt/sammas_storage/bergland-lab/Daphnia_HWE/hwe_stat.Rdata")
  hwe.stat[,py:=paste(pond, year, sep=".")]

### load precomputed file (HWE_simulations.prep.rivanna.R makes this file)
  load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
  sc[,year:=tstrsplit(clone, "_")[[2]]]
  sc[,pond:=gsub("Dcat", "DCat", pond)]

### load in snp.set *getTargetSNPs.analysis.workstation.R makes this file)
  load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/snpSet.Rdata")

### open genotype file
  genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### function to return table
  getGeno <- function(py.list.snps   = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
                      py.list.clones = py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019)),
                      sc.test=c("M", "OO.93"),
                      snps = snp.set) {

    ### define pond to use
      #py.list.snps <- list(c("D8"), c(2016, 2017, 2018, 2019))
      #py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019))
      #sc.test=c("M", "OO.93")
      #py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019))
      #snps <- snp.set


    ### define py.i for SNPs
      py.i.snps <-paste(paste(py.list.snps[[1]], collapse="."), paste(py.list.snps[[2]], collapse="."), sep=".")
      py.i.clones <-paste(paste(py.list.clones[[1]], collapse="."), paste(py.list.clones[[2]], collapse="."), sep=".")

      print(paste(py.i.snps, py.i.clones, sep="  .  "))

    ### first subsample clones to use in assay
      clones <- sc[Species=="Pulex"][year%in%py.list.clones[[2]]][pond%in%py.list.clones[[1]], list(clone=sample(clone, size=1)), list(sc.uniq)]
      clones[,pond:=tstrsplit(clone, "_")[[3]]]
      clones[,year:=tstrsplit(clone, "_")[[2]]]

    ### SNPs to use
      #sites <- snp.set[py==py.i.snps][(zw)]
      sites <- snp.set[py==py.i.snps]

    ### which of our clones arefocal clones
      setkey(clones, clone)
      setkey(sc, clone)
      AB <- merge(clones, sc)[sc.uniq.x%in%sc.test]

    ### get genotype data
      seqResetFilter(genofile)
      seqSetFilter(genofile, variant.id=sites$variant.id, sample.id=AB$clone)

      genomat <- as.data.table(t(seqGetData(genofile, "$dosage")))
      setnames(genomat, names(genomat), AB[match(seqGetData(genofile, "sample.id"), AB$clone)]$sc.uniq.x)

      genomat[,variant.id:=seqGetData(genofile, "variant.id")]
      genomat[,chr:=seqGetData(genofile, "chromosome")]
      genomat[,pos:=seqGetData(genofile, "position")]

    ### wide to long
      gl <- melt(genomat, id.vars=c("variant.id", "chr", "pos"))
      setnames(gl, c("variable", "value"), c("sc", "geno"))

      gl[,sc := factor(sc, sc.test)]

    ### tack in block info
      setkey(gl, variant.id)
      setkey(sites, variant.id)

      glo <- merge(gl, sites)

    ### return
      glo[,py.clones:=py.i.clones]
      glo[,py.snps:=py.i.snps]
      return(glo)
  }

### tab geno function
  tabGeno <- function(tmp) {
    ### iterate through block ids
      tmp.ag <- foreach(i=unique(tmp$block.id), .errorhandling="remove", .combine="rbind")%dopar%{
          foreach(zw.i=c(TRUE, FALSE), .errorhandling="remove", .combine="rbind")%do%{
            tt <- tmp[block.id==i][zw==zw.i]

            tab <- table(tt[geno!=0]$sc, tt[geno!=0]$geno)

            #fet <- fisher.test(tab)

            data.table(block.id=i, n.SNPs=length(unique(tt$variant.id)),
                       KB=tt$KB[1], chr=tt$chr.x[1], BP1=tt$BP1[1], BP2=tt$BP2[1],
                       A.0=sum(tt[sc==levels(sc)[1]]$geno==0, na.rm=T),
                       A.1=sum(tt[sc==levels(sc)[1]]$geno==1, na.rm=T),
                       A.2=sum(tt[sc==levels(sc)[1]]$geno==2, na.rm=T),
                       B.0=sum(tt[sc==levels(sc)[2]]$geno==0, na.rm=T),
                       B.1=sum(tt[sc==levels(sc)[2]]$geno==1, na.rm=T),
                       B.2=sum(tt[sc==levels(sc)[2]]$geno==2, na.rm=T),
                       zw=zw.i)

          }
      }
      tmp.ag[,freqA_het:=A.1/(A.0 + A.1 + A.2)]
      tmp.ag[,freqB_het:=B.1/(B.0 + B.1 + B.2)]

    ### return
      return(tmp.ag)
  }

### plotting function
  plotGeno <- function(tmp.ag) {


    ### some formatting
      tmp.ag[,congruent:= (freqA_het>.50 & freqB_het==0)]
      #table(tmp.ag$congruent)
      #summary(glm(congruent~n, tmp.ag, family=binomial()))
      #summary(glm(congruent~KB, tmp.ag, family=binomial()))
      #summary(glm(congruent~chr, tmp.ag, family=binomial()))

      tmp.ag.ag <- tmp.ag[,list(n=length(KB)), list(freqB_het, freqA_het, zw, x)]

      g1 <- ggplot(data=tmp.ag.ag, aes(x=freqA_het, y=freqB_het, size=n,  color=as.factor(x))) +
      geom_point() + xlim(0,1) + ylim(0,1) + facet_grid(~zw)

      return(g1)
  }

### D8
    D8.sc <- unique(sc[pond=="D8"]$sc.uniq)[c(1:2)]

    D8.o <- foreach(i=1:(length(D8.sc)-1), .errorhandling="remove")%dopar%{
      o <- foreach(j=(i+1):length(D8.sec), .errorhandling="remove")%do%{
        geno <- getGeno(py.list.snps   = py.list.snps <- list(c("D8"), c(2016, 2017, 2018, 2019)),
                              py.list.clones = py.list.clones <- list(c("D8"), c(2016, 2017, 2018, 2019)),
                              sc.test=c(D8.sc[i], D8.sc[j]))
        geno[,scA:=D8.sc[i]]
        geno[,scB:=D8.sc[j]]
        geno
      }
      rbindlist(o)
    }
    D8.o <- rbindlist(D8.o)
    save(D8.o, file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/D8.pairwise.Rdata")

    





    D8.D8.tab <- tabGeno(D8.D8.geno)
    D8.D8.tab[block.id%in%D8.D8.tab[zw==T & freqA_het==1 & freqB_het==0]$block.id, x:=0]
    D8.D8.tab[is.na(x), x:=1]

    plotGeno(D8.D8.tab[KB>10][n.SNPs>10]) + geom_hline(yintercept=0.05)

    D8.D8.tab[zw==T & freqA_het==1 & freqB_het==0, x:="a"]
    D8.D8.tab[zw==F & freqB_het<.05, x:="b"]


    w <- dcast(D8.D8.tab, block.id~zw, value.var="x")

    x1 <- D8.D8.tab[zw==T][freqA_het==1 & freqB_het==0]$block.id
    x2 <- D8.D8.tab[zw==F][freqB_het<.05]$block.id

    D8.D8.tab[block.id%in%intersect(x1, x2), foo:=T]
    D8.D8.tab[is.na(foo), foo:=F]

    ggplot(data=D8.D8.tab[KB>10][n.SNPs>10], aes(x=foo, y=log10(n.SNPs/KB), group=foo, color=as.factor(zw))) +
    geom_point()



        ggplot(data=D8.D8.tab, aes(x=foo, y=log10(KB), group=foo)) +
        geom_point()
  ### DBunk on D8:A/B
    DBunk.D8.geno <- getGeno(py.list.snps   = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
                    py.list.clones = py.list.clones <- list(c("D8"), c(2016, 2017, 2018, 2019)),
                    sc.test=c("A", "B"))


    DBunk.D8.tab <- tabGeno(DBunk.D8.geno)
    p1 <- plotGeno(DBunk.D8.tab[zw==T][n>10])
    p2 <- plotGeno(DBunk.D8.tab[zw==F][n>10])
    plot_grid(p1, p2)

  ### D8 on DCat:A*/B*
    D8.DCat.geno <- getGeno(py.list.snps   = py.list.snps <- list(c("D8"), c(2016, 2017, 2018, 2019)),
                    py.list.clones = py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019)),
                    sc.test=c("M", "OO.93"))


    D8.DCat.tab <- tabGeno(D8.DCat.geno)
    p1 <- plotGeno(D8.DCat.tab[zw==T])
    p2 <- plotGeno(D8.DCat.tab[zw==F])
    plot_grid(p1, p2)

    x1 <- D8.DCat.tab[zw==T][freqA_het==1 & freqB_het==0]$block.id
    x2 <- D8.DCat.tab[zw==F][freqB_het<.01]$block.id

    intersect(x1, x2)

  ### DBunk on DCat:A*/B*
    DBunk.DCat.geno <- getGeno(py.list.snps   = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
                    py.list.clones = py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019)),
                    sc.test=c("M", "OO.93"))


    DBunk.DCat.tab <- tabGeno(DBunk.DCat.geno)
    p1 <- plotGeno(DBunk.DCat.tab[zw==T])
    p2 <- plotGeno(DBunk.DCat.tab[zw==F])
    plot_grid(p1, p2)


table(c(D8.DCat.tab[freqA_het==1 & freqB_het==0][KB>10]$block.id,
        DBunk.DCat.tab[freqA_het==1 & freqB_het==0][KB>10]$block.id))


table(D8.DCat.tab$freqA_het==1 & D8.DCat.tab$freqB_het==0)
table(DBunk.DCat.tab$freqA_het==1 & DBunk.DCat.tab$freqB_het==0)










      plot(freqA_het~freqB_het, tmp.ag[n>20])




    ### DBunk on D8
      tmp <- getGeno(py.list.snps   = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
              py.list.clones = py.list.clones <- list(c("D8"), c(2016, 2017, 2018, 2019)),
              sc.test=c("A", "B"))
      fisher.test(table(tmp[geno!=0]$sc, tmp[geno!=0]$geno))

      tmp.ag <- foreach(i=unique(tmp$block.id), .errorhandling="remove", .combine="rbind")%dopar%{
        tab <- table(tmp[block.id==i][geno!=0]$sc, tmp[block.id==i][geno!=0]$geno)

        fet <- fisher.test(tab)

        data.table(block.id=i, or=fet$estimate, p=fet$p.value, n=sum(tab))
      }

### DCat clones
  ### D8 on DCat

    tmp <- getGeno(py.list.snps   = list(c("D8"), c(2016, 2017, 2018, 2019)),
            py.list.clones = list(c("DCat"), c(2016, 2017, 2018, 2019)),
            sc.test=c("M", "OO.93"))
    fisher.test(table(tmp[geno!=0]$sc, tmp[geno!=0]$geno))

    tmp.ag <- foreach(i=unique(tmp$block.id), .errorhandling="remove", .combine="rbind")%dopar%{

      tt <- tmp[block.id==i]

      tab <- table(tt[geno!=0]$sc, tt[geno!=0]$geno)

      fet <- fisher.test(tab)

      data.table(block.id=i, or=fet$estimate, p=fet$p.value, n=sum(tab), KB=tt$KB[1], BP1=tt$BP1[1], BP2=tt$BP2[1],
                 A.0=sum(tt[sc==levels(sc)[1]]$geno==0, na.rm=T),
                 A.1=sum(tt[sc==levels(sc)[1]]$geno==1, na.rm=T),
                 A.2=sum(tt[sc==levels(sc)[1]]$geno==2, na.rm=T),
                 B.0=sum(tt[sc==levels(sc)[2]]$geno==0, na.rm=T),
                 B.1=sum(tt[sc==levels(sc)[2]]$geno==1, na.rm=T),
                 B.2=sum(tt[sc==levels(sc)[2]]$geno==2, na.rm=T))

    }
    tmp.ag[,freqA_het:=A.1/(A.0 + A.1 + A.2)]
    tmp.ag[,freqB_het:=B.1/(B.0 + B.1 + B.2)]


    plot(freqA_het~freqB_het, tmp.ag[n>20])













    tmp.ag <- tmp[NSNPS>10,list(or=fisher.test(table(sc, geno))$estimate), list(block.id)]


    ### simple test. A is heterozygote, B is homozygote.
      table(glo[block.id==1]$sc, glo[block.id==1]$geno)













    ### second, calculate Euclidian distance from point of highest difference based on the simulations
      maxDiff <- hwe.ag.m.ag[py==py.i][which.max(diff.mu)]
      hwe.stat[py==py.i, euclid.dist := (fAA - maxDiff$fAA)^2 + (fAa - maxDiff$fAa)^2 + (faa - maxDiff$faa)^2]

    ### extract out SNPs
      sites <- hwe.stat[py==py.i & euclid.dist<=1e-4]

    ### which of our clones are A &B
      setkey(clones, clone)
      setkey(sc, clone)
      AB <- merge(clones, sc)[sc.uniq.x%in%sc.test]

    ### get genotype data
      seqResetFilter(genofile)
      seqSetFilter(genofile, variant.id=sites$variant.id, sample.id=AB$clone)

      genomat <- as.data.table(t(seqGetData(genofile, "$dosage")))
      setnames(genomat, names(genomat), AB[match(seqGetData(genofile, "sample.id"), AB$clone)]$sc.uniq.x)

      genomat[,variant.id:=seqGetData(genofile, "variant.id")]
      genomat[,chr:=seqGetData(genofile, "chromosome")]
      genomat[,pos:=seqGetData(genofile, "position")]

    ### wide to long
      gl <- melt(genomat, id.vars=c("variant.id", "chr", "pos"))
      setnames(gl, c("variable", "value"), c("sc", "geno"))

      gl[,sc := factor(sc, sc.test)]

    ### simple test. A is heterozygote, B is homozygote.
      table(gl[geno>0]$sc, gl[geno>0]$geno)
  }

### run function
  ### D8 clones
    genoTab(py.list.snps   = py.list.snps <- list(c("D8"), c(2016, 2017, 2018, 2019)),
            py.list.clones = py.list.clones <- list(c("D8"), c(2016, 2017, 2018, 2019)),
            sc.test=c("A", "B"))

    genoTab(py.list.snps   = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
            py.list.clones = py.list.clones <- list(c("D8"), c(2016, 2017, 2018, 2019)),
            sc.test=c("A", "B"))

  ### Bunk clones
    genoTab(py.list.snps   = py.list.snps <- list(c("D8"), c(2016, 2017, 2018, 2019)),
            py.list.clones = py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019)),
            sc.test=c("M", "OO.93"))

    genoTab(py.list.snps   = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
            py.list.clones = py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019)),
            sc.test=c("M", "OO.93"))
