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

### LD clumping

  ld.blocks <- foreach(py.list=list(list(c("D8"), c(2016, 2017, 2018, 2019)),
                                    list(c("DBunk"), c(2016, 2017, 2018, 2019))))%do%{
    ### define py.i
      py.i <-paste(paste(py.list[[1]], collapse="."), paste(py.list[[2]], collapse="."), sep=".")
      print(py.i)

    ### first subsample clones
      clones <- subsampClone(sc.dt=sc[Species=="Pulex"][year%in%py.list[[2]]], use.pond=py.list[[1]])

    ### second, get out allele frequencies of subsampled clones
     #  seqSetFilter(genofile,
     #               sample.id=clones$clone,
     #               variant.id=snp.dt[(final.use)]$id)

     #  clones.af <- data.table(id=seqGetData(genofile, "variant.id"),
     #                          af=seqAlleleFreq(genofile, .progress=T, parallel=F))

     #  setkey(clones.af, id)
     #  setkey(snp.dt, id)

     #  clones.af <- merge(clones.af, snp.dt)[af<0 & af<1] #[af<1/length(seqGetData(genofile, "sample.id")) &
     #                                                     #af<1-1/length(seqGetData(genofile, "sample.id"))]

   ## third, calculate Euclidian distance from point of highest difference based on the simulations
      maxDiff <- hwe.ag.m.ag[py==py.i][which.max(diff.mu)]
      hwe.stat[py==py.i, euclid.dist := (fAA - maxDiff$fAA)^2 + (fAa - maxDiff$fAa)^2 + (faa - maxDiff$faa)^2]

   ### fourth, write make fake bigsnpr file
      seqResetFilter(genofile)
      seqSetFilter(genofile, variant.id=hwe.stat[euclid.dist<1e-4]$variant.id, sample.id=clones$clone)

      seqGDS2VCF(genofile, vcf.fn="/tmp/tmp.vcf", info.var=character(0), fmt.var=character(0))

    ### SNP clumping
      system(paste("plink1.9 --double-id --vcf /tmp/tmp.vcf --allow-extra-chr --blocks no-pheno-req --blocks-max-kb 5000"))

    ### load in data
      blocks <- fread("~/plink.blocks.det")
      blocks[,py:=py.i]

    ### return
      blocks

  }
  ld.blocks <- rbindlist(ld.blocks)

### save
  save(ld.blocks, file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/ldBlocks.Rdata")

### a little bit of summary statistics. Wohooo!
  ld.blocks[,list(chr=CHR[which.max(KB)], start=BP1[which.max(KB)], stop=BP2[which.max(KB)], n=NSNPS[which.max(KB)]), list(py)]
  ld.blocks[,list(chr=CHR[which.max(NSNPS)], start=BP1[which.max(NSNPS)], stop=BP2[which.max(NSNPS)], n=NSNPS[which.max(NSNPS)]), list(py)]



### test of overlaping ranges
  A <- makeGRangesFromDataFrame(data.frame(chr=ld.blocks[py=="D8.2016.2017.2018.2019"][KB>0]$CHR,
                                      start=ld.blocks[py=="D8.2016.2017.2018.2019"][KB>0]$BP1,
                                      end=ld.blocks[py=="D8.2016.2017.2018.2019"][KB>0]$BP2),
                           start.field="start", end.field="end")

  B <- makeGRangesFromDataFrame(data.frame(chr=ld.blocks[py=="DBunk.2016.2017.2018.2019"][KB>0]$CHR,
                                     start=ld.blocks[py=="DBunk.2016.2017.2018.2019"][KB>0]$BP1,
                                     end=ld.blocks[py=="DBunk.2016.2017.2018.2019"][KB>0]$BP2),
                          start.field="start", end.field="end")

  genome <- getGenomeAndMask(genome=snp.dt[(final.use), list(start=min(pos), end=max(pos)), list(chr)], mask=NA)

  pt <- overlapPermTest(A=A, B=B, genome=genome$genome, ntimes=1000)
  plot(pt)

### make plot
  ol <- as.data.table(overlapRegions(A=A, B=B,  min.bases=50000))
  setnames(ol, "chr", "seqnames")

  ggplot() +
  geom_segment(data=as.data.table(genome$genome),
                aes(x=start, xend=end, y=1, yend=1)) +
  facet_wrap(~seqnames, scales="free_x") +
  geom_rect(data=ol,
            aes(xmin=startA, xmax=endA, ymin=-1, ymax=4), fill = "grey") +
  geom_segment(data=as.data.table(A),
                aes(x=start, xend=end, y=2, yend=2), size=2, color="red") +
  geom_segment(data=as.data.table(B),
                aes(x=start, xend=end, y=3, yend=3), size=2, color="blue")
