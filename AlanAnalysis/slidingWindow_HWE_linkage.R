### libraries
  library(SeqArray)
  library(SNPRelate)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(SeqVarTools)
  library(ggtern)

### functions
  ### funciton to radomly subset one clone per superclone
    subsampClone <- function(sc.dt, n=1, use.pond="DBunk") {
      sc.samp <- sc.dt[,list(clone=sample(clone, size=n)), list(sc.uniq)]
      sc.samp[,pond:=tstrsplit(clone, "_")[[3]]]
      sc.samp[,year:=tstrsplit(clone, "_")[[2]]]
      return(sc.samp[pond%in%use.pond])
    }


### make SeqArray object
  #seqVCF2GDS(vcf.fn="/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.vcf",
  #           out.fn="/mnt/spicy_3/AlanDaphnia/vcf/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds")

### load data
  ### filtered SNPs
    #load("/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/finalsnpstousewSimoids_20190430.Rdata")
    load("/mnt/spicy_3/Karen/201620172018FinalMapping/snpsetfilteredformissing_20190423.Rdata")
    use <- data.table(id=snpsetfilteredformissing, use=T)
    setkey(use, id)

  ### superclones
    sc <- fread("/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/Superclones20161718withlowcoverageindupdated_20190501")
    sc[,pond := tstrsplit(clone, "_")[[3]]]
    sc[,sc.uniq := SC]
    sc[SC=="OO", sc.uniq:=paste(SC, SCnum, sep=".")]

  ### open GDS object
    #genofile <- snpgdsOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds")
    #genofile <- snpgdsOpen("/mnt/ssd/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds")
    #genofile <- snpgdsOpen("/mnt/ssd/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds")
    genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

  ### make SNP table
    ### import and merge with filtering file
      snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                           pos=seqGetData(genofile, "position"),
                           id=seqGetData(genofile, "variant.id"))
      setkey(snp.dt, id)
      snp.dt <- merge(snp.dt, use, all.x=T, all.y=T)

    ### subset on the good chromosomes
      snp.dt.ag <- snp.dt[,list(nSNPs=length(pos)), list(chr)]
      snp.dt.ag[nSNPs>3000, use.chr:=T]
      snp.dt.ag[is.na(use.chr), use.chr:=F]

      setkey(snp.dt.ag, chr)
      setkey(snp.dt, chr)

      snp.dt <- merge(snp.dt, snp.dt.ag, all.x=T)
      snp.dt[is.na(use), use:=F]

    ### final SNP filter
      snp.dt[,final.use := use & use.chr]

### get allele frequencies for DBunk, downsampled to one per superclone
  clones <- subsampClone(sc.dt=sc, use.pond="DBunk")

  seqSetFilter(genofile,
               sample.id=clones[year==2017]$clone,
               variant.id=snp.dt[(final.use)]$id)

  clones.af <- data.table(id=seqGetData(genofile, "variant.id"),
                          af=seqAlleleFreq(genofile, .progress=T, parallel=F))

### save
  save(clones, clones.af, snp.dt, sc, file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")

### load
  load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
  genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds", allow.duplicate=TRUE)

  setkey(clones.af, id)
  setkey(snp.dt, id)

  clones.af <- merge(clones.af, snp.dt)[af>.01 & af<.99]

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







### generate HWE & F across the genome
  seqResetFilter(genofile)
  seqSetFilter(genofile, sample.id=clones$clone, variant.id=clones.af$id)

  hwe <- as.data.table(hwe(genofile))
  setkey(hwe, variant.id)

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
