# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


### libraries
  library(data.table)
  #library(ggplot2)
  #library(patchwork)
  library(SeqArray)
  library(SeqVarTools)


setwd("/scratch/aob2x/daphnia_hwe_sims")
### load F1 phenotype data here made here: `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")

### load superclone file to get SC.uniq
  sc <- fread("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  sc <- sc[Nonindependent==0 ]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]
  sc[,pond:=toupper(population)]


  ### libraries
    library(data.table)
    library(foreach)

  ### convert to GDS
    #vcf.fn <- "/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf"
    #gds.fn <- "/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.gds"
    #seqVCF2GDS(vcf.fn, gds.fn)

    gds.fn <- "/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds"

  ### open genofile
    genofile <- seqOpen(gds.fn)
    snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                         pos=seqGetData(genofile, "position"),
                         id=seqGetData(genofile, "variant.id"),
                         numAlleles=seqNumAllele(genofile),
                         ref=seqGetData(genofile, "$ref"),
                         alt=seqGetData(genofile, "$alt"))
    setkey(snp.dt, chr, pos)

  ### clone metadata file
    sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

    sc <- sc[Nonindependent==0 ]

    sc[,SC.uniq:=SC]
    sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]
    sc[,pond:=toupper(population)]

    ### sc per year
      sc.peryear <- sc[,list(.N), list(year, SC.uniq, pond)]

    ### hard filtering of SC
      sc.ag <- sc[LabGenerated==F & Species=="pulex", list(clone=clone[which.max(medrd)][1]), list(SC.uniq, Species)]
      #sc.ag[,pond:=toupper(pond)]

  ### load poolseq data
    load("/project/berglandlab/alan/gprime_peaks.replicates.250K.05.Rdata")
    setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))
    setkey(gprime, chr, pos,rep)

  ### function
    polarize <- function(i, window) {
      #window<- 10000; i<-10
      chr.i<-peaks[i]$CHROM; start<-peaks[i]$posMaxGprime-window; stop<-peaks[i]$posMaxGprime+window; rep.i=peaks[i]$rep

      pool.tmp <- gprime[J(data.table(chr=chr.i, pos=start:stop, rep=rep.i, key="chr,pos,rep")), nomatch=0]
      setkey(pool.tmp, chr, pos)
      setkey(snp.dt, chr, pos)
      unique(snp.dt[pool.tmp]$id)
      seqSetFilter(genofile, variant.id=unique(snp.dt[pool.tmp]$id))

      tmp <- as.data.table(getGenotype(genofile))
      tmp[,sample.id:=seqGetData(genofile, "sample.id")]

      dat.phase <- melt(tmp, id.vars="sample.id", variable.name="variant.id", value.name="geno")
      dat.phase[,allele1:=tstrsplit(geno, "\\|")[[1]]]
      dat.phase[,allele2:=tstrsplit(geno, "\\|")[[2]]]

      # dat.phase[sample.id=="April_2017_D8_213"] ## A
      # dat.phase[sample.id=="April_2017_D8_151"] ## C

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
