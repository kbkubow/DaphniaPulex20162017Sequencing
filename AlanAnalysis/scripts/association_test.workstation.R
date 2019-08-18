### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(ggtern)
  library(viridis)
  library(bigstatsr)
  library(SeqArray)
  library(bigsnpr)
  library(regioneR)
  library(doMC)
  registerDoMC(20)

  ### load phenotype data
    pheno <- fread("/mnt/spicy_3/AlanDaphnia/male_female_pheno/male_female_dorset.csv")

    pheno[,clone.id:=clone_ID]
    pheno[,clone.id:=gsub("AD", "D", clone.id)]
    pheno[,clone.id:=gsub("AW", "W", clone.id)]

    pheno[,clone:=unlist(sapply(pheno$clone.id, function(x)  sc[grepl(x, clone)]$clone))]


### load HWE summary data (HWE_simulations.gather.rivanna.R makes this file)
  load(file="/mnt/sammas_storage/bergland-lab/Daphnia_HWE/hwe.ag.m.ag")

### open HWE & genotype dist statistic data
  ### cp /scratch/aob2x/daphnia_hwe_sims/hwe_stat.Rdata /nv/vol186/bergland-lab/Daphnia_HWE/hwe_stat.Rdata
  load(file="/mnt/sammas_storage/bergland-lab/Daphnia_HWE/hwe_stat.Rdata")
  hwe.stat[,py:=paste(pond, year, sep=".")]

### load precomputed file (HWE_simulations.prep.rivanna.R makes this file)
  load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
  sc[,year:=tstrsplit(clone, "_")[[2]]]

### open genotype file
  genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### downsample function
  subsampClone <- function(sc.dt, n=1, use.pond="DBunk") {
    sc.samp <- sc.dt[,list(clone=sample(clone, size=n)), list(sc.uniq)]
    sc.samp[,pond:=tstrsplit(clone, "_")[[3]]]
    sc.samp[,year:=tstrsplit(clone, "_")[[2]]]
    return(sc.samp[pond%in%use.pond])
  }

### first subsample clones
  clones <- subsampClone(sc.dt=sc[clone%in%unique(pheno$clone)], use.pond=c("D8", "DBunk"))

  setkey(pheno, clone)
  setkey(clones, clone)
  sub.pheno <- merge(pheno, clones)

  sub.pheno.ag <- sub.pheno[,list(nMale=mean(total_males, na.rm=T),
                                  fracMale=round(mean(total_males))/round(mean(total_individuals)),
                                  n=round(mean(total_individuals)),
                                  epp=total_ephippia>0),
                              list(clone, pond=pond.x)]

### second, calculate Euclidian distance from point of highest difference based on the simulations
  sites <- foreach(py.i = c("D8.2016.2017.2018.2019", "DBunk.2017"), .combine="rbind")%do%{

    maxDiff <- hwe.ag.m.ag[py==py.i][which.max(diff.mu)]
    hwe.stat[py==py.i, euclid.dist := (fAA - maxDiff$fAA)^2 + (fAa - maxDiff$fAa)^2 + (faa - maxDiff$faa)^2]

  ### extract out SNPs
    hwe.stat[py==py.i & euclid.dist<=1e-4]
  }
  setkey(sites, chr, pos)
  sites <- sites[!duplicated(sites)]

### get genotype data
  seqSetFilter(genofile, variant.id=sites$variant.id, sample.id=sub.pheno.ag$clone)

  genomat <- as.data.table(t(seqGetData(genofile, "$dosage")))
  setnames(genomat, names(genomat), sub.pheno.ag[match(seqGetData(genofile, "sample.id"), sub.pheno.ag$clone)]$clone)

  genomat[,variant.id:=seqGetData(genofile, "variant.id")]
  genomat[,chr:=seqGetData(genofile, "chromosome")]
  genomat[,pos:=seqGetData(genofile, "position")]

### wide to long
  gl <- melt(genomat, id.vars=c("variant.id", "chr", "pos"))
  setnames(gl, c("variable", "value"), c("clone", "geno"))


### merge with phenotype data forsimple test
  gl.ag <- gl[geno!=0,list(mu=mean(geno, na.rm=T)), list(clone)]

  setkey(gl.ag, clone)
  setkey(sub.pheno.ag, clone)

  m <- merge(gl.ag, sub.pheno.ag)

  summary(lm(fracMale~mu, m))

### GWASy type thing
  setkey(gl, clone)
  mgl <- merge(gl, sub.pheno.ag, allow.cartesian=T)

  setkey(mgl, variant.id)
  o <- foreach(i=unique(mgl$variant.id), .errorhandling="remove", .combine="rbind")%dopar%{
    ### i <- 1547618
    print(i)
    tmp <- mgl[J(i)]

      #t1 <- lm(nMale~geno, tmp[geno!=0][pond%in%c("D8", "DBunk")])
      t1 <- glm(fracMale~geno, tmp[geno!=0][pond%in%c("D8", "DBunk")],
            weights=tmp[geno!=0][pond%in%c("D8", "DBunk")]$n, family="binomial")

    data.table(variant.id=i, i=which(i==unique(mgl$variant.id)), p=summary(t1)$coef[2,4])
  }

  plot(I(-log10(p))~i, o)
