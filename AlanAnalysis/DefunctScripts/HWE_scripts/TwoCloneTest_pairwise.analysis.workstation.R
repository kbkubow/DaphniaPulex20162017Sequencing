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

            #tab <- table(tt[geno!=0]$sc, tt[geno!=0]$geno)

            #fet <- fisher.test(tab)

            data.table(block.id=i, n.SNPs=length(unique(tt$variant.id)),
                       KB=tt$KB[1], chr=tt$chr.x[1], BP1=tt$BP1[1], BP2=tt$BP2[1],
                       A.0=sum(tt[sc==unique(sc)[1]]$geno==0, na.rm=T),
                       A.1=sum(tt[sc==unique(sc)[1]]$geno==1, na.rm=T),
                       A.2=sum(tt[sc==unique(sc)[1]]$geno==2, na.rm=T),
                       B.0=sum(tt[sc==unique(sc)[2]]$geno==0, na.rm=T),
                       B.1=sum(tt[sc==unique(sc)[2]]$geno==1, na.rm=T),
                       B.2=sum(tt[sc==unique(sc)[2]]$geno==2, na.rm=T),
                       zw=zw.i, scA=tt$scA[1], scB=tt$scB[1])

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

      tmp.ag.ag <- tmp.ag[,list(n=length(KB)), list(freqB_het, freqA_het, zw)]

      g1 <- ggplot(data=tmp.ag.ag, aes(x=freqA_het, y=freqB_het, size=n)) +
      geom_point() + xlim(0,1) + ylim(0,1) + facet_grid(~zw)

      return(g1)
  }

### D8
  D8.sc <- unique(sc[pond=="D8"]$sc.uniq)
  D8.o <- foreach(i=1:(length(D8.sc)-1), .errorhandling="remove")%dopar%{
    o <- foreach(j=(i+1):length(D8.sc), .errorhandling="remove")%do%{
      print(paste(i, j, sep=" . "))
      geno <- getGeno(py.list.snps = py.list.snps <- list(c("D8"), c(2016, 2017, 2018, 2019)),
                            py.list.clones = py.list.clones <- list(c("D8"), c(2016, 2017, 2018, 2019)),
                            sc.test=c(D8.sc[i], D8.sc[j]))
      geno[,scA:=D8.sc[i]]
      geno[,scB:=D8.sc[j]]
      write.csv(geno,  file=paste("/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/D8_", i, "_", j, ".csv", sep=""))
      return("foo")
    }
    o
  }

### DBunk
  DBunk.sc <- unique(sc[pond=="DBunk"]$sc.uniq)
  DBunk.o <- foreach(i=1:(length(DBunk.sc)-1), .errorhandling="remove")%dopar%{
    o <- foreach(j=(i+1):length(DBunk.sc), .errorhandling="remove")%do%{
      print(paste(i, j, sep=" . "))
      geno <- getGeno(py.list.snps = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
                            py.list.clones = py.list.clones <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
                            sc.test=c(DBunk.sc[i], DBunk.sc[j]))
      geno[,scA:=DBunk.sc[i]]
      geno[,scB:=DBunk.sc[j]]
      write.csv(geno,  file=paste("/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/DBunk_", i, "_", j, ".csv", sep=""))
      return("foo")

    }
    (o)
  }


### load in for post-processing
### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(ggplot2)
  library(viridis)

### aggregate
  fl <- system("ls /mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/D*_*_*.csv", intern=T)

  D8.o <- foreach(i=fl[grepl("D8", fl)])%dopar%{
    print(paste(which(i==fl[grepl("D8", fl)]), " / ", length(fl[grepl("D8", fl)])))
    #i <- "/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/D8_10_11.csv"
    dat <- fread(i)

    dat[,sc.f := LETTERS[as.numeric(as.factor(sc))]]

    dat.tab <- dat[,list(n.SNPs=length(variant.id), KB=mean(KB, na.rm=T), chr=chr.x[1], BP1=BP1[1], BP2=BP2[1],
              n.GT=c(sum(geno==0, na.rm=T), sum(geno==1, na.rm=T), sum(geno==2, na.rm=T)),
              GT=c(0,1,2)),
         list(block.id, zw, sc, sc.f)]

    dat.tab.ag <- dat.tab[,list(n=sum(!is.na(KB)), KB=mean(KB, na.rm=T)), list(block.id)]

    dat.tab <- dat.tab[block.id%in%dat.tab.ag[n==12]$block.id]
    dat.tab[,pond:=tstrsplit(last(tstrsplit(i, "/")), "_")[[1]]]
    dat.tab[,pair.i:=which(i==fl[grepl("D8", fl)])]

    #plotGeno(dat.tab)

    dat.tab
  }
  D8.o <- rbindlist(D8.o)
  save(D8.o, file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/D8_ag.Rdata")


  DBunk.o <- foreach(i=fl[grepl("DBunk", fl)])%dopar%{
    print(paste(which(i==fl[grepl("DBunk", fl)]), " / ", length(fl[grepl("DBunk", fl)])))
    #i <- "/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/D8_10_11.csv"
    dat <- fread(i)

    dat[,sc.f := LETTERS[as.numeric(as.factor(sc))]]

    dat.tab <- dat[,list(n.SNPs=length(variant.id), KB=mean(KB, na.rm=T), chr=chr.x[1], BP1=BP1[1], BP2=BP2[1],
              n.GT=c(sum(geno==0, na.rm=T), sum(geno==1, na.rm=T), sum(geno==2, na.rm=T)),
              GT=c(0,1,2)),
         list(block.id, zw, sc, sc.f)]

    dat.tab.ag <- dat.tab[,list(n=sum(!is.na(KB)), KB=mean(KB, na.rm=T)), list(block.id)]

    dat.tab <- dat.tab[block.id%in%dat.tab.ag[n==12]$block.id]
    dat.tab[,pond:=tstrsplit(last(tstrsplit(i, "/")), "_")[[1]]]
    dat.tab[,pair.i:=which(i==fl[grepl("DBunk", fl)])]

    #plotGeno(dat.tab)

    dat.tab
  }
  DBunk.o <- rbindlist(DBunk.o)
  save(DBunk.o, file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/DBunk_ag.Rdata")


### process aggregates
  load("D8_ag.Rdata")
  #D8.o <- DBunk.o
  D8.o.ag <- D8.o[KB>20,list(het.zw=n.GT[GT==1 & (zw)]/sum(n.GT[(zw)], na.rm=T),
                             het.nzw=n.GT[GT==1 & !(zw)]/sum(n.GT[!(zw)], na.rm=T),
                             het=n.GT[GT==1]/sum(n.GT, na.rm=T),
                             n=sum(n.GT)),
                        list(sc, pair.i, block.id)]

  D8.o.ag.ag <- D8.o.ag[,list(het.zw=mean(het.zw, na.rm=T),
                              het.nzw=mean(het.nzw, na.rm=T),
                              het=mean(het, na.rm=T),
                              n=mean(n)),
                         list(sc, block.id)]

plot(het.zw~het, D8.o.ag.ag[sc=="A"])
plot(het.zw~het, D8.o.ag.ag[sc=="B"])



  ### for sorting on SC
    D8.o.ag.ag.ag <- D8.o.ag.ag[,list(n=sum(het.mu>.95, na.rm=T) - sum(het.mu<.05, na.rm=T),
                                      meanHet=mean(het.mu, na.rm=T), meanHom=mean(het.mu<.05, na.rm=T)),
                                  list(sc, zw)]


    D8.o.ag.ag[,sc:=factor(sc, levels=D8.o.ag.ag.ag[(zw)][order(meanHet)]$sc)]

  ### for sorting on block
    D8.o.ag.ag.ag <- D8.o.ag.ag[,list(n=sum(het.mu>.95, na.rm=T),
                                    meanHet=mean(het.mu, na.rm=T), meanHom=mean(het.mu<.05, na.rm=T)),
                                list(block.id, zw)]
    D8.o.ag.ag[,block.id:=factor(block.id, levels=D8.o.ag.ag.ag[(zw)][order(meanHet)]$block.id)]

  ggplot(data=D8.o.ag.ag, aes(x=block.id, y=sc, fill=het.mu)) + geom_tile() + facet_wrap(~zw)


  load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/D8_ag.Rdata")

  o.ag <- o[KB>1,list(fracHet=n.GT[GT==1]/sum(n.GT, na.rm=T), KB=mean(KB), nSNPs=sum(n.GT)), list(block.id, zw, sc, pair.i, sc.f)]

  o.ag.w <- dcast(o.ag, value.var=c("fracHet", "KB"), block.id  + pair.i + zw ~ sc.f)
  o.ag.w.ag <- o.ag.w[,list(n=length(KB_A)),
                       list(pair.i, zw, fracHet_A, fracHet_B)]



  ggplot(data=o.ag.w.ag, aes(x=fracHet_A, y=fracHet_B, size=n)) + geom_hex() + facet_grid(zw~.)



#### quick association test on blocks

  ### prepare pheno data
    ### load SC data

      load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
      sc[,year:=tstrsplit(clone, "_")[[2]]]


    ### load phenotype data
      pheno <- fread("/mnt/spicy_3/AlanDaphnia/male_female_pheno/male_female_dorset.csv")

      pheno[,clone.id:=clone_ID]
      pheno[,clone.id:=gsub("AD", "D", clone.id)]
      pheno[,clone.id:=gsub("AW", "W", clone.id)]

      pheno[,clone:=unlist(sapply(pheno$clone.id, function(x)  sc[grepl(x, clone)]$clone))]


      pheno <- merge(pheno, sc, by="clone")


  load("/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/D8_ag.Rdata")
  D8.o.ag <- D8.o[KB>20,list(het=n.GT[GT==1]/sum(n.GT, na.rm=T)), list(sc, pair.i, block.id, zw)]
  D8.o.ag.ag <- D8.o.ag[,list(het.mu=mean(het, na.rm=T)), list(sc, block.id, zw)]
  m <- merge(D8.o.ag.ag, pheno, by.x="sc", by.y="SC", allow.cartesian=T)
  m.ag <- m[,list(frac=mean(total_males/total_individuals, na.rm=T), het=mean(het.mu)),
             list(sc, block.id, zw)]


  load("/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/DBunk_ag.Rdata")
  DBunk.o.ag <- DBunk.o[KB>20,list(het=n.GT[GT==1]/sum(n.GT, na.rm=T)), list(sc, pair.i, block.id, zw)]
  DBunk.o.ag <- DBunk.o[KB>20,list(het=n.GT[GT==1]/sum(n.GT, na.rm=T)), list(sc, pair.i, block.id, zw)]
  DBunk.o.ag.ag <- DBunk.o.ag[,list(het.mu=mean(het, na.rm=T)), list(sc, block.id, zw)]
  m <- merge(DBunk.o.ag.ag, pheno, by.x="sc", by.y="SC", allow.cartesian=T)
  m.ag <- m[,list(frac=mean(total_males/total_individuals, na.rm=T), het=mean(het.mu)),
             list(sc, block.id, zw)]


  m.ag.ag <- m.ag[,list(cor=cor(frac, het, use="complete")), list(zw, block.id)]






















  o.ag.w <- dcast(o.ag[zw==F], value.var=c("fracHet"), pair.i ~ zw + sc.f + block.id )
  levelplot(t(as.matrix(na.omit(o.ag.w))[,-1]))

  mat <-
  hc <- hclust()


  o.ag <- o[KB>20,list(fracHet_A=n.GT[GT==1 & sc.f=="A"]/sum(n.GT[sc.f=="A"], na.rm=T),
                      fracHet_B=n.GT[GT==1 & sc.f=="B"]/sum(n.GT[sc.f=="B"], na.rm=T),
                      KB=mean(KB), nSNPs=sum(n.GT),
                      scA=sc[sc.f=="A"][1], scB=sc[sc.f=="B"][1]),
                 list(block.id, zw, pair.i)]

 o.ag[,delta := (fracHet_A - fracHet_B)]


 o.ag.ag <- o.ag[,list(median_delta=median(delta, na.rm=T), KB=mean(KB)), list(block.id, zw)]

 tmp <- o.ag[,list(ratio=nSNPs[(zw)]/nSNPs[(!zw)], KB=mean(KB)), list(block.id, pair.i)]


 fuck <- foreach(i=unique(o.ag$block.id), .errorhandling="remove")%do%{
   print(i)
   i<-7
   tmp <- o.ag[zw==T][block.id==i]
   tmp.ag <- tmp[,list(n=length(KB)), list(fracHet_A, fracHet_B)]
   ggplot(data=tmp.ag, aes(x=fracHet_A, y=fracHet_B, size=n)) + geom_point()

   setkey(tmp, scA, scB)
   mat <- matrix(NA, nrow=length(unique(tmp$scA))+1, ncol=length(unique(tmp$scA))+1)
   mat[lower.tri(mat, diag=F)] <- tmp$delta

   mat <- Matrix::forceSymmetric(mat,uplo="L")

   diag(mat) <- 0

   #heatmap.2(as.matrix((mat)))


   #mc <- cmdscale(as.matrix(abs(mat)), 1)
   as.matrix(mat)
}





 plot(mc[,1]~mc[,2])

 ggplot(data=o.ag[zw==T][block.id==1][order(delta)], aes(x=scA, y=scB, fill=abs(delta))) + geom_tile() + scale_color_viridis() + scale_fill_viridis()

   o.ag.w <- dcast(o.ag[zw==F], value.var=c("delta"), pair.i ~ block.id )
   levelplot(as.matrix(o.ag.w)[,-1])
