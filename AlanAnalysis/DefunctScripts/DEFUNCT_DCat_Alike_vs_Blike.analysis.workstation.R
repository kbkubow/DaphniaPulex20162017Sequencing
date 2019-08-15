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
  sc[,pond:=gsub("Dcat", "DCat", pond)]


### open genotype file
  genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### function to return table
  genoTab <- function(py.list.snps   = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
                      py.list.clones = py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019)),
                      sc.test=c("M", "OO.93")) {

    ### define pond to use
      #py.list.snps <- list(c("D8"), c(2016, 2017, 2018, 2019))
      #py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019))

      #py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019))

    ### define py.i for SNPs
      py.i <-paste(paste(py.list.snps[[1]], collapse="."), paste(py.list.snps[[2]], collapse="."), sep=".")
      print(py.i)

    ### first subsample clones to use in assay
      clones <- sc[Species=="Pulex"][year%in%py.list.clones[[2]]][pond%in%py.list.clones[[1]], list(clone=sample(clone, size=1)), list(sc.uniq)]
      clones[,pond:=tstrsplit(clone, "_")[[3]]]
      clones[,year:=tstrsplit(clone, "_")[[2]]]

    ### second, calculate Euclidian distance from point of highest difference based on the simulations
      maxDiff <- hwe.ag.m.ag[py==py.i][which.max(diff.mu)]
      hwe.stat[py==py.i, euclid.dist := (fAA - maxDiff$fAA)^2 + (fAa - maxDiff$fAa)^2 + (faa - maxDiff$faa)^2]

    ### extract out SNPs
      sites <- hwe.stat[py==py.i & euclid.dist<=1e-4]

    ### which of our clones are A &B
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

    ### simple test. A is heterozygote, B is homozygote.
      table(gl[geno>0]$sc, gl[geno>0]$geno)
  }

### run function
  ### D8 clones
    genoTab(py.list.snps   = py.list.snps <- list(c("D8"), c(2016, 2017, 2018, 2019)),
            py.list.clones = py.list.clones <- list(c("D8"), c(2016, 2017, 2018, 2019)),
            sc.test=c("A", "B"))

    genoTab(py.list.snps   = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
            py.list.clones = py.list.clones <- list(c("D8"), c(2016, 2017, 2018, 2019)),
            sc.test=c("A", "B"))

  ### Bunk clones
    genoTab(py.list.snps   = py.list.snps <- list(c("D8"), c(2016, 2017, 2018, 2019)),
            py.list.clones = py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019)),
            sc.test=c("M", "OO.93"))

    genoTab(py.list.snps   = py.list.snps <- list(c("DBunk"), c(2016, 2017, 2018, 2019)),
            py.list.clones = py.list.clones <- list(c("DCat"), c(2016, 2017, 2018, 2019)),
            sc.test=c("M", "OO.93"))



1233/500
1434/410
