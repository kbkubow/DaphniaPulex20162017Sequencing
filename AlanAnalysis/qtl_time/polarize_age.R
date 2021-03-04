### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  
### load poolseq data
  load("/project/berglandlab/alan/gprime_peaks.replicates.250K.05.Rdata")
  setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))
  setkey(gprime, chr, pos,rep)

### convert to GDS
  vcf.fn <- "/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf"
  gds.fn <- "/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.gds"

  #seqVCF2GDS(vcf.fn, gds.fn)

### open genofile
  genofile <- seqOpen(gds.fn)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       ref=seqGetData(genofile, "$ref"),
                       alt=seqGetData(genofile, "$alt"))
  setkey(snp.dt, chr, pos)


### for each qtl
  qtl_age <- foreach(qtl.n=1:14)%do%{
    #qtl.n <- 10
    message(qtl.n)
    st <- sign(gprime[chr==peaks$CHROM[qtl.n] & pos==peaks$posMaxGprime[qtl.n] & rep==peaks$rep[qtl.n]]$deltaSNP)

    if(st==1) male_limited_allele=0; male_plus_allele=1
    if(st==-1) male_limited_allele=1; male_plus_allele=0

    ### load correct Tca file for male limited allele
      fn.male_limited <- paste("/project/berglandlab/jcbnunez/Daphnia_Age_Est/TCA_figure/State_",
                  male_limited_allele, "/",
                  "QTL_", qtl.n, ".recode_mcmc_list.RDATA", sep="")

      fn.male_plus <- paste("/project/berglandlab/jcbnunez/Daphnia_Age_Est/TCA_figure/State_",
                  male_plus_allele, "/",
                  "QTL_", qtl.n, ".recode_mcmc_list.RDATA", sep="")

      load(fn.male_limited)
      male_limited_age=median(mcmc.output$t.chain[,1])

      load(fn.male_plus)
      male_male_plus_age=median(mcmc.output$t.chain[,1])

    ### get derived state
      seqSetFilter(genofile, variant.id=snp.dt[chr==peaks$CHROM[qtl.n] & pos==peaks$posMaxGprime[qtl.n]]$id,
                  sample.id=c("pulicaria", "obtusa"))
      anc <- seqGetData(genofile, "$dosage_alt")

      data.table(qtl=qtl.n,
                  male_limited_age=male_limited_age,
                  male_plus_age=male_male_plus_age,
                  ancestral=ifelse(male_plus_allele==anc[1], "male_plus", "male_limited"),
                  male_limited_allele=male_limited_allele,
                  pulicaria=anc[1,], obtusa=anc[2,])
  }
  qtl_age <- rbindlist(qtl_age)
