#ijob -c1 -p standard -A berglandlab
# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### set working directory
setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

#
### load meta-data file
  samps <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### load GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

### load in filter file
  snpFilter <- fread("snpsvarpulexpresentinhalf_table_20200623")

### make snp.dt
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       key="chr")
  setkey(snpFilter, chr, pos)
  setkey(snp.dt, chr, pos)

  snp.dt <- merge(snpFilter, snp.dt)



### 1/2N PNPS function
  pnps.fun <- function(SC.i="A", MAC=1, sample.n=50) {
    #SC.i="A"; MAC=1; sample.n=10
    seqResetFilter(genofile)
    seqSetFilter(genofile,
                 sample.id=sample(samps[population%in%c("D8")][SC==SC.i][Nonindependent==0]$clone, sample.n),
                 variant.id=snpFilter$variant.ids)

    message("getting allele counts")
    snp.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                         alleleCount=seqAlleleCount(genofile, ref.allele=1L))



    ### get annotation
    message("getting annotations")

      snp.dt <- snp.dt[alleleCount==MAC]
      seqSetFilter(genofile, variant.id=snp.dt$variant.id)

      tmp <- seqGetData(genofile, "annotation/info/ANN")
      len1 <- tmp$length
      len2 <- tmp$data

      snp.dt1 <- data.table(len=rep(len1, times=len1),
                            ann=len2,
                            id=rep(snp.dt$variant.id, times=len1))

    # Extracting data between the 2nd and third | symbol
      snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]

    # Collapsing additional annotations to original SNP vector length
      snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=",")),
                            list(variant.id=id)]

      snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]

    # merge

      m <- merge(snp.dt, snp.dt1.an, by="variant.id")

    ### get ref and alt allele frequencies
      alleleDepth <- seqGetData(genofile, "annotation/format/AD")
      refDepth <-

    ### summarize
      m.ag <- m[alleleCount==MAC][,list(pn=sum(col=="missense_variant"), ps=sum(col=="synonymous_variant"), SC=SC.i), list(MAC=alleleCount)]

    ### return
      return(m.ag)
  }


  ### 1/2N PNPS function
    refFreq.fun <- function(SC.i="A", MAC=1, sample.n=50) {
      #SC.i="OO"; MAC=1; sample.n=10
      seqResetFilter(genofile)
      seqSetFilter(genofile,
                   sample.id=sample(samps[population%in%c("D8")][SC==SC.i][Nonindependent==0]$clone, sample.n),
                   variant.id=snpFilter$variant.ids)

      message("getting allele counts")
      snp.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                           alleleCount=seqAlleleCount(genofile, ref.allele=1L))


      ### get annotation
      message("getting annotations")

      snp.dt <- snp.dt[alleleCount==MAC]
      seqSetFilter(genofile, variant.id=snp.dt$variant.id)
s
      ### summarize
        m.ag <- m[alleleCount==MAC][,list(pn=sum(col=="missense_variant"), ps=sum(col=="synonymous_variant"), SC=SC.i), list(MAC=alleleCount)]

      ### return
        return(m.ag)
    }


### run
  pnps.out <- foreach(j=1:50, .combine="rbind")%do%{
    o <- foreach(i=c("A", "C", "OO"), .combine="rbind")%do%{
      pnps.fun(SC.i=i, sample.n=25)
    }
    o[,boot:=j]
  }

### SAVE
  save(pnps.out, file="/nv/vol186/bergland-lab/alan/pnps_zurich.Rdata")

### download
  # scp aob2x@rivanna.hpc.virginia.edu:/nv/vol186/bergland-lab/alan/pnps_zurich.Rdata .
  scp aob2x@rivanna.hpc.virginia.edu:~/within_A.Rdata ~/.

### load, plot, done.
  library(ggplot2)
  library(data.table)
  library(cowplot)

  load("pnps_zurich.Rdata")
  load("~/within_A.Rdata")
  pnps.out[,pnps:=pn/ps]
  pnps.out.ag <- pnps.out[,list(pnps=mean(pnps), lci=quantile(pnps, .025), uci=quantile(pnps, .975)), list(SC)]

  ggplot(data=pnps.out.ag) +
  geom_errorbar(aes(x=SC, ymin=lci, ymax=uci), width=.1) +
  geom_point(aes(x=SC, y=pnps)) +
  scale_x_discrete(labels=c("A" = "A", "C" = "B",
                              "OO" = "Singleton \n clones")) +
  ylab("pN/pS") + xlab("") + ggtitle("1/2N mutations")
