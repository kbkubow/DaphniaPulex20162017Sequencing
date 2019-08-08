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

### get SNP list (uses same distance parameter as LD_clump script)

  ### define pond to use
    py.list<- list(c("D8"), c(2016, 2017, 2018, 2019))

  ### define py.i
    py.i <-paste(paste(py.list[[1]], collapse="."), paste(py.list[[2]], collapse="."), sep=".")
    print(py.i)

  ### first subsample clones
    clones <- subsampClone(sc.dt=sc[Species=="Pulex"][year%in%py.list[[2]]], use.pond=py.list[[1]])

  ### second, calculate Euclidian distance from point of highest difference based on the simulations
    maxDiff <- hwe.ag.m.ag[py==py.i][which.max(diff.mu)]
    hwe.stat[py==py.i, euclid.dist := (fAA - maxDiff$fAA)^2 + (fAa - maxDiff$fAa)^2 + (faa - maxDiff$faa)^2]

  ### extract out SNPs
    sites <- hwe.stat[py==py.i & euclid.dist<=1e-4]

  ### which of our clones are A &B
    setkey(clones, clone)
    setkey(sc, clone)
    AB <- merge(clones, sc)[SC%in%c("A", "B")]

  ### get genotype data
    seqSetFilter(genofile, variant.id=sites$variant.id, sample.id=AB$clone)

    genomat <- as.data.table(t(seqGetData(genofile, "$dosage")))
    setnames(genomat, names(genomat), AB[match(seqGetData(genofile, "sample.id"), AB$clone)]$SC)

    genomat[,variant.id:=seqGetData(genofile, "variant.id")]
    genomat[,chr:=seqGetData(genofile, "chromosome")]
    genomat[,pos:=seqGetData(genofile, "position")]

  ### wide to long
    gl <- melt(genomat, id.vars=c("variant.id", "chr", "pos"))
    setnames(gl, c("variable", "value"), c("sc", "geno"))

  ### simple test. A is heterozygote, B is homozygote.
    table(gl[geno>0]$sc, gl[geno>0]$geno)
