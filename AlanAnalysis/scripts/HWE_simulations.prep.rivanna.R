### libraries
  library(SeqArray)
  #library(SNPRelate)
  library(data.table)
  library(foreach)
  #library(doMC)
  registerDoMC(20)
  library(SeqVarTools)
  #library(ggtern)

### functions
  ### funciton to radomly subset one clone per superclone
    subsampClone <- function(sc.dt, n=1, use.pond="DBunk") {
      sc.samp <- sc.dt[,list(clone=sample(clone, size=n)), list(sc.uniq)]
      sc.samp[,pond:=tstrsplit(clone, "_")[[3]]]
      sc.samp[,year:=tstrsplit(clone, "_")[[2]]]
      return(sc.samp[pond%in%use.pond])
    }

### make SeqArray object
  #seqVCF2GDS(vcf.fn="/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.vcf",
  #           out.fn="/mnt/spicy_3/AlanDaphnia/vcf/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds")

### load data
  ### filtered SNPs

    #load("/mnt/spicy_3/Karen/201620172018FinalMapping/snpsetfilteredformissing_20190423.Rdata") ### this goes with A
    #use <- data.table(id=snpsetfilteredformissing, use=T)

    load("/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/finalsnpstousewSimoids_20190430.Rdata") ### this goes with version D; use this
    use <- data.table(id=finalsnpstousewSimoids, use=T)

    setkey(use, id)

  ### superclones
    #sc <- fread("/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/Superclones20161718withlowcoverageindupdated_20190501")
    #sc <- fread("/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/Superclones20161718withlowcoverageindupdated_20190802.csv")
    #sc <- fread("/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/Superclones20161718updated_20190802.csv")
    sc <- fread("/mnt/spicy_3/AlanDaphnia/vcf/Superclones20161718updated_20190802.csv")
    sc[,pond := tstrsplit(clone, "_")[[3]]]
    sc[,sc.uniq := SC]
    sc[SC=="OO", sc.uniq:=paste(SC, SCnum, sep=".")]

  ### open GDS object
    #genofile <- snpgdsOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds")
    #genofile <- snpgdsOpen("/mnt/ssd/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds")
    #genofile <- snpgdsOpen("/mnt/ssd/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.gds")
    #genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
    #genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
    genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

  ### make SNP table
    ### import and merge with filtering file
      seqSetFilter(genofile, sample.id=sc[Species=="Pulex"]$clone)
      snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                           pos=seqGetData(genofile, "position"),
                           id=seqGetData(genofile, "variant.id"))
      setkey(snp.dt, id)
      snp.dt <- merge(snp.dt, use, all.x=T, all.y=T)

    ### subset on the good chromosomes
      snp.dt.ag <- snp.dt[,list(nSNPs=length(pos)), list(chr)]
      snp.dt.ag[nSNPs>3000, use.chr:=T]
      snp.dt.ag[is.na(use.chr), use.chr:=F]

      setkey(snp.dt.ag, chr)
      setkey(snp.dt, chr)

      snp.dt <- merge(snp.dt, snp.dt.ag, all.x=T)
      snp.dt[is.na(use), use:=F]

    ### final SNP filter
      snp.dt[,final.use := use & use.chr]

### save
  save(snp.dt, sc, file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")


######## pick up here #######
### gets ported to Rivanna under file name: HWE_simulations_rivanna.R


  ### libraries
    library(SeqArray)
    library(data.table)
    library(foreach)
    #library(doMC)
    #registerDoMC(20)
    library(SeqVarTools)
    #library(ggtern)
    #library(viridis)

  ### pull in SLURM ARRAY BATCH ID
    	args.vec <- as.numeric(commandArgs(trailing=T)) + 1
    #    args.vec <- 1


  ### funciton to radomly subset one clone per superclone
    subsampClone <- function(sc.dt, n=1, use.pond="DBunk") {
      sc.samp <- sc.dt[,list(clone=sample(clone, size=n)), list(sc.uniq)]
      sc.samp[,pond:=tstrsplit(clone, "_")[[3]]]
      sc.samp[,year:=tstrsplit(clone, "_")[[2]]]
      return(sc.samp[pond%in%use.pond])
    }

  ### load precomuted file
    #load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
    load(file="/scratch/aob2x/daphnia_hwe_sims/subFiles.Rdata")
    sc[,year:=tstrsplit(clone, "_")[[2]]]

  ### open GDS object
    #genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds", allow.duplicate=TRUE)
    #genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
    genofile <- seqOpen("/scratch/aob2x/daphnia_hwe_sims/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

  ### downsampled to one per superclone, get allele frequencies
    table(sc$pond, sc$year, sc$Species)
    #unique(sc[Species=="Pulex"][!grepl("W", pond)]$pond)

    clones.l <-  foreach(pond=list("D8", "DBunk", c("D8", "DBunk")))%dopar%{
        foreach(year.i=list(2017, c(2016, 2017, 2018, 2019)))%dopar%{

          print(paste(paste(pond, collapse="."), paste(year.i, collapse="."), sep=" / "))

          clones <- subsampClone(sc.dt=sc[Species=="Pulex"][year%in%year.i], use.pond=pond)

          seqSetFilter(genofile,
                       sample.id=clones$clone,
                       variant.id=snp.dt[(final.use)]$id)

          clones.af <- data.table(id=seqGetData(genofile, "variant.id"),
                                  af=seqAlleleFreq(genofile, .progress=T, parallel=F))

          setkey(clones.af, id)
          setkey(snp.dt, id)

          clones.af <- merge(clones.af, snp.dt)[af>0 & af<1] #[af>1/length(seqGetData(genofile, "sample.id")) &
                                                             #af<1-1/length(seqGetData(genofile, "sample.id"))]


          o <-  list()
          o$pond <- pond
          o$year <-  year.i
          o$clones <- clones
          o$clones.af  <- clones.af

          return(o)
        }
    }
    #save(clones.l, file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/clones_l.D8.DBunk.Rdata")

  ### generate HWE & F across the genome
    hwe.stat <- foreach(pond.i=1:length(clones.l),  .combine="rbind")%dopar%{
      foreach(year.i=1:length(clones.l[[pond.i]]),  .combine="rbind")%dopar%{
        print(paste(paste(clones.l[[pond.i]][[year.i]]$pond, collapse="."),
                    paste(clones.l[[pond.i]][[year.i]]$year, collapse="."),
                    sep=" / "))

        seqResetFilter(genofile)
        seqSetFilter(genofile,
                     sample.id=clones.l[[pond.i]][[year.i]]$clones$clone,
                     variant.id=clones.l[[pond.i]][[year.i]]$clones.af$id)

        hwe <- as.data.table(hwe(genofile))
        setkey(hwe, variant.id)

        hwe[,fAA:=nAA/(nAA+nAa+naa)]
        hwe[,fAa:=nAa/(nAA+nAa+naa)]
        hwe[,faa:=naa/(nAA+nAa+naa)]

        hwe <- merge(hwe, snp.dt, by.x="variant.id", by.y="id")
        hwe[,pond:=paste(clones.l[[pond.i]][[year.i]]$pond, collapse=".")]
        hwe[,year:=paste(clones.l[[pond.i]][[year.i]]$year, collapse=".")]

        return(hwe)
      }
    }

    save(hwe.stat, file="/scratch/aob2x/daphnia_hwe_sims/hwe_stat.Rdata")
