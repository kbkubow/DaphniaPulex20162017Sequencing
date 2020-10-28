# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(foreach)
  library(data.table)


### load data
  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/lmer4qtl/", "AxC", full.name=T)
  o <- foreach(fn.i=fn)%do%{
    #fn.i <- fn[1]
    message(fn.i)
    load(fn.i)
    return(lmer.gwas)
  }
  o <- rbindlist(o)

### summarize
  o.ag <- o[,list(minp=min(p.z, na.rm=T), maxc=max(chisq)), list(term, perm)]

  o.ag[perm>0, list(q=quantile(minp, .05)), list(term)]



### libraries
  library(SeqArray)

  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile))
  setkey(snp.dt, id)

  snp.dt[id==976200]
  seqSetFilter(genofile, variant.id=976200)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")


  target <- data.table(dosage=seqGetData(genofile, "$dosage"),
                      clone=seqGetData(genofile, "sample.id"))

  target <- merge(target, sc, "clone")

  table(na.omit(target[Species=="pulex"][population=="D8"][,list(dosage=mean(dosage.V1, na.rm=T)), list(SCnum)])$dosage)


)



ggplot() +
geom_line(data=gprime, aes(x=POS, y=G, color=CHROM)) +
geom_point(data=gprime[POS==3968943], aes(x=POS, y=G)) +
facet_grid(~CHROM)
