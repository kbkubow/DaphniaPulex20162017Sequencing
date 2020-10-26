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









ggplot() +
geom_line(data=gprime, aes(x=POS, y=G, color=CHROM)) +
geom_point(data=gprime[POS==3968943], aes(x=POS, y=G)) +
facet_grid(~CHROM)
