#ijob -c1 -p standard -A berglandlab
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### load genotyping data
  #load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
  load(file="/scratch/aob2x/daphnia_hwe_sims/subFiles.Rdata")
  sc[,year:=tstrsplit(clone, "_")[[2]]]

  ### open GDS object
  #genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds", allow.duplicate=TRUE)
  #genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
  genofile <- seqOpen("/scratch/aob2x/daphnia_hwe_sims/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

### get differences between A&B
  genodat <- foreach(sc.i=c("A", "B"))%do%{
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sc[SC==sc.i][pond=="D8"]$clone)

    tmp <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     variant.id=seqGetData(genofile, "variant.id"),

                     sc=sc.i,
                     af= 1 - seqAlleleFreq(genofile))
  }
  genodat <- rbindlist(genodat)

  genodat.w <- dcast(genodat, chr + pos + variant.id ~ sc)


### pooled data
  load("/nv/vol186/bergland-lab/alan/totalADRDlongall.Rdata")

  geno[,effRD:=floor(RRD)]
  geno[,effPA:=round(propalt*effRD)/effRD]

  geno.w <- dcast(geno[pond=="D8"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

### merge
  setkey(genodat.w, chr, pos)
  setkey(geno.w, chr, pos)

  m <- merge(genodat.w, geno.w, all.x=T, all.y=T)

### are pooled D8 samples basically F1s?
  m <- na.omit(m)
  m[,m.hat:=(effPA_D8Male1 * effRD_D8Male1 + effPA_D8Male2 * effRD_D8Male2) / (effRD_D8Male1 + effRD_D8Male2)]
  m[,f.hat:=(effPA_D8PE1 * effRD_D8PE1 + effPA_D8PE2 * effRD_D8PE2) / (effRD_D8PE1 + effRD_D8PE2)]
  m[!is.na(A),A.geno := unlist(sapply(m[!is.na(A)]$A, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[!is.na(B),B.geno := unlist(sapply(m[!is.na(B)]$B, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]



  m.ag <- m[,list(mu=c(mean(m.hat), mean(f.hat)),
          sd=c(sd(m.hat), sd(f.hat)),
          pool=c("m", "f")),
      list(A.geno, B.geno)]
  m.ag[,exp.fq:=(A.geno+B.geno)/2]


  save(m.ag, file="/nv/vol186/bergland-lab/alan/mf_expectation.Rdata")



  library(ggplot2)
  library(data.table)

  load("/mnt/sammas_storage/bergland-lab/alan/mf_expectation.Rdata")

  F1.plot <- ggplot(data=m.ag, aes(x=pool, y=mu)) +
  geom_hline(aes(yintercept=exp.fq/2), linetype="solid", color="red", size=.75) +
  geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd), width=.2, size=.75) +
  geom_point(size=4, fill="white", shape=21, stroke=1) +
  facet_grid(~A.geno+B.geno) +
  ylab("Frequencuy") +
  xlab("Pool") +
  scale_x_discrete(labels=c("m" = "Male", "f" = "Female")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
