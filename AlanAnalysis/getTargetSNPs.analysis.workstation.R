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

### open HWE & genotype dist statistic data
  ### cp /scratch/aob2x/daphnia_hwe_sims/hwe_stat.Rdata /nv/vol186/bergland-lab/Daphnia_HWE/hwe_stat.Rdata
  load(file="/mnt/sammas_storage/bergland-lab/Daphnia_HWE/hwe_stat.Rdata")
  hwe.stat[,py:=paste(pond, year, sep=".")]

### load HWE summary data (HWE_simulations.gather.rivanna.R makes this file)
  load(file="/mnt/sammas_storage/bergland-lab/Daphnia_HWE/hwe.ag.m.ag")

### load in LD clump data
  load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/ldBlocks.Rdata")

### load precomputed file (HWE_simulations.prep.rivanna.R makes this file)
  load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
  sc[,year:=tstrsplit(clone, "_")[[2]]]

### open genotype file
  genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### make LD clump data long to intersect with snp.dt
  ld.blocks[py=="D8.2016.2017.2018.2019", block.id:=rank(-KB, ties="first")]
  ld.blocks[py=="DBunk.2016.2017.2018.2019", block.id:=rank(-KB, ties="first")]

  ld.blocks.long <- ld.blocks[,list(chr=CHR, pos=c(BP1:BP2), KB=KB, NSNPS=NSNPS, BP1=BP1, BP2=BP2), list(py, block.id)]

### intersect
  setkey(ld.blocks.long, chr, pos)
  setkey(snp.dt, chr, pos)

  m <- merge(snp.dt, ld.blocks.long)
  setnames(m, "id", "variant.id")

### get SNP sets
  subsampClone <- function(sc.dt, n=1, use.pond="DBunk") {
    sc.samp <- sc.dt[,list(clone=sample(clone, size=n)), list(sc.uniq)]
    sc.samp[,pond:=tstrsplit(clone, "_")[[3]]]
    sc.samp[,year:=tstrsplit(clone, "_")[[2]]]
    return(sc.samp[pond%in%use.pond])
  }

  snp.set <- foreach(py.list=list(list(c("D8"), c(2016, 2017, 2018, 2019)),
                                    list(c("DBunk"), c(2016, 2017, 2018, 2019))))%do%{
    ### define py.i
      py.i <-paste(paste(py.list[[1]], collapse="."), paste(py.list[[2]], collapse="."), sep=".")
      print(py.i)

    ### first subsample clones
      clones <- subsampClone(sc.dt=sc[Species=="Pulex"][year%in%py.list[[2]]], use.pond=py.list[[1]])

   ## third, calculate Euclidian distance from point of highest difference based on the simulations
      maxDiff <- hwe.ag.m.ag[py==py.i][which.max(diff.mu)]
      hwe.stat[py==py.i, euclid.dist := (fAA - maxDiff$fAA)^2 + (fAa - maxDiff$fAa)^2 + (faa - maxDiff$faa)^2]

   ### get output basic intersection sets
      zw.snps <- data.table(py=py.i,
                 variant.id=hwe.stat[py==py.i & euclid.dist<1e-4]$variant.id,
                 zw=T)

    ### merge
      setkey(m, variant.id, py)
      setkey(zw.snps, variant.id, py)

      m.tmp <- merge(m, zw.snps, all.x=T)[py==py.i]
      m.tmp[is.na(zw), zw:=F]

    ### get frequencies for SNPs
      setkey(hwe.stat, variant.id, py)
      setkey(m.tmp, variant.id, py)
      m.out <- merge(m.tmp, hwe.stat[,c("variant.id", "py", "nAA", "nAa", "naa", "afreq", "f"), with=F])

    ### return
      return(m.out)
  }
  snp.set <- rbindlist(snp.set)

### save
  save(snp.set, file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/snpSet.Rdata")
