# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


### libraries
  library(data.table)
  #library(ggplot2)
  #library(patchwork)
  library(SeqArray)
  library(SeqVarTools)
  library(foreach)

  setwd("/scratch/aob2x/daphnia_hwe_sims")

### load F1 phenotype data here made here: `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")
  male.ag <- male[!is.na(gr),list(propmale=sum(Males)/sum(NewTotal),
                                  N=sum(NewTotal)),
                    list(clone, gr)]
  male.ag[,se:=sqrt((propmale*(1-propmale))/N)]

### load superclone file to get SC.uniq
  sc <- fread("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  sc <- sc[Nonindependent==0 ]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]
  sc[,pond:=toupper(population)]

### load genotype data (from DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.convert_output.vitPath.R)
  pio <- fread(file= "/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all.vitPath.pio.csv")
  pio[,phase.fold.geno:=phase.geno]
  pio[phase.fold.geno%in%c(12, 21), phase.fold.geno:=12]

### trim to samples used
  setkey(pio, clone)
  pio <- pio[J(unique(c(male$clone)))]

### load poolseq data
  load("/project/berglandlab/alan/gprime_peaks.replicates.250K.05.Rdata")
  setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))
  setkey(gprime, chr, pos,rep)

### function
  polarize <- function(i, window) {
    #window<- 10000; i<-10
    message(i)
    chr.i<-peaks[i]$CHROM; start<-peaks[i]$posMaxGprime-window; stop<-peaks[i]$posMaxGprime+window; rep.i=peaks[i]$rep

    pool.tmp <- gprime[J(data.table(chr=chr.i, pos=start:stop, rep=rep.i, key="chr,pos,rep")), nomatch=0]
    setkey(pool.tmp, chr, pos)
    setkey(pio, chr, pos)
    unique(pio[pool.tmp]$id)

    dat.phase <- merge(pool.tmp, pio)[,c("chr", "pos", "clone", "deltaSNP", "phase.fold.geno"), with=F]
    dat.phase[,allele1:=substr( phase.fold.geno, 0,1)]
    dat.phase[,allele2:=substr( phase.fold.geno, 2,3)]
    dat.phase

    dat.phase <- melt(dat.phase[,-"phase.fold.geno",with=F],
                    id.vars=c("clone", "chr", "pos", "deltaSNP"),
                    value.vars=c("allele1", "allele2"),
                    value.name="allele")

    dat.phase[sign(deltaSNP)== -1 &  allele==2, concord:="male"]
    dat.phase[sign(deltaSNP)== -1 &  allele==1, concord:="pe"]
    dat.phase[sign(deltaSNP)==  1 &  allele==2, concord:="pe"]
    dat.phase[sign(deltaSNP)==  1 &  allele==1, concord:="male"]

    m.ag <- dcast(dat.phase, chr+pos+clone~variable, value.var="concord")
    m.ag.ag <- m.ag[,list(n.male_male=sum(allele1=="male" & allele2=="male")/length(allele1),
                  n.male_pe=sum((allele1=="male" & allele2=="pe") | (allele2=="male" & allele1=="pe"))/length(allele1),
                  n.pe_pe=sum(allele1=="pe" & allele2=="pe")/length(allele1)),
              list(clone)]

    m.ag.ag[,qtl:=i]
    m.ag.ag

  }

### iterage and summarize
  qtl.polar.f1 <- foreach(i=1:14, .combine="rbind")%do%polarize(i, window=15000)
  qtl.polar.f1.ag <- qtl.polar.f1[,list(geno=c("male_male", "male_pe", "pe_pe")[which.max(c(n.male_male, n.male_pe, n.pe_pe))]), list(clone, qtl)]

### merge
  f1.pool.merge <- merge(male.ag, qtl.polar.f1.ag, by="clone")


### save
  save(f1.pool.merge, file="DaphniaPulex20162017Sequencing/AlanFigures/Figure4/f1_pool_polar.Rdata")

### load
  library(data.table)
  library(ggplot)

  load()


























    setkey(qtl.polar, SC.uniq)
    setkey(sc.peryear, SC.uniq)
    qtl.polar <- merge(qtl.polar, sc.peryear)

              save(qtl.polar, file="~/qtl_polar.Rdata")













      dat.phase <- melt(dat.phase[,c("sample.id", "variant.id", "allele1", "allele2")],
                      id.vars=c("sample.id", "variant.id"),
                      value.vars=c("allele1", "allele2"),
                      value.name="allele")

      dat.phase[,variable:=gsub("allele", "", variable)]
      dat.phase[,haplo:=paste(sample.id, variable, sep=".")]
      setnames(dat.phase, "variant.id", "id")
      dat.phase[,id:=as.numeric(as.character(id))]

      setkey(dat.phase, id)
      setkey(snp.dt, id)


      dat.phase <- merge(dat.phase, snp.dt)
      dat.phase[,rep:=rep.i]

      setkey(dat.phase, chr, pos, rep)
      setkey(gprime, chr, pos, rep)

      m <- merge(dat.phase, gprime)

      m[sign(deltaSNP)== -1 &  allele==0, concord:="male"]
      m[sign(deltaSNP)== -1 &  allele==1, concord:="pe"]
      m[sign(deltaSNP)==  1 &  allele==0, concord:="pe"]
      m[sign(deltaSNP)==  1 &  allele==1, concord:="male"]
      m <- merge(m, sc.ag, by.x="sample.id", by.y="clone")

      # dcast(m[SC.uniq=="C"], chr+pos+SC.uniq~variable, value.var="concord")
      # dcast(m[SC.uniq=="A"], chr+pos~variable, value.var="concord")

      m.ag <- dcast(m, chr+pos+SC.uniq~variable, value.var="concord")
      setnames(m.ag, c("1","2"), c("allele1", "allele2"))
      m.ag.ag <- m.ag[,list(n.male_male=sum(allele1=="male" & allele2=="male")/length(allele1),
                    n.male_pe=sum((allele1=="male" & allele2=="pe") | (allele2=="male" & allele1=="pe"))/length(allele1),
                    n.pe_pe=sum(allele1=="pe" & allele2=="pe")/length(allele1)),
                list(SC.uniq)]
      setkey(m.ag.ag, SC.uniq)
      setkey(sc.ag, SC.uniq)

      m.ag.ag <- merge(m.ag.ag, sc.ag)
      m.ag.ag[,pond:=toupper(tstrsplit(clone, "_")[[3]])]
      m.ag.ag[,qtl:=i]
      m.ag.ag

    }


    qtl.polar <- foreach(i=1:14, .combine="rbind")%do%polarize(i, window=15000)

    setkey(qtl.polar, SC.uniq)
    setkey(sc.peryear, SC.uniq)
    qtl.polar <- merge(qtl.polar, sc.peryear)

    save(qtl.polar, file="~/qtl_polar.Rdata")






























### tack SC.uniq onto the male F1 data
  male <- merge(male, sc[,c("clone", "SC.uniq"), with=F], by="clone")


  tmp <- male[OneLiterPheno==1][,list(propmale=mean(propmale)), list(SC.uniq)]

  setkey(qtl.polar, SC.uniq)
  setkey(tmp, SC.uniq)

  qtl.polar.sc <- qtl.polar[J(tmp), nomatch=0]
  qtl.polar.sc[,male.freq:=2*n.male_male + n.male_pe]
  qtl.polar.sc.l <- melt(qtl.polar.sc[,c("SC.uniq", "qtl", "n.male_male", "n.male_pe", "n.pe_pe", "propmale")], value.var=c("qtl", "SC.uniq", "propmale"), measure.var=c("n.male_male", "n.male_pe", "n.pe_pe"))
  qtl.polar.sc.l.ag <- qtl.polar.sc.l[,list(geno=variable[which.max(value)], propmale=mean(propmale)), list(SC.uniq, qtl)]

  summ <- qtl.polar.sc.l.ag[,list(.N, propmale=mean(propmale)), list(geno, SC.uniq)]

plot(propmale~N, summ[geno=="n.male_male"])
  ggplot(data=qtl.polar.sc) +
  geom_line(aes(y=male.freq, x=SC.uniq, group=qtl), position=position_dodge(width=0.3), color="grey", size=.5) +
  geom_point(aes(y=male.freq, x=SC.uniq, color=SC.uniq, group=qtl), position=position_dodge(width=0.3))
