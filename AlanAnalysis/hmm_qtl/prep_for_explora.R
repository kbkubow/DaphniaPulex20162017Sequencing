### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)
  library(ggplot2)

### pooled data
  load("/nv/vol186/bergland-lab/alan/totalADRDlongall.Rdata")

  geno[,effRD:=floor(RRD)]
  geno[,effPA:=round(propalt*effRD)/effRD]


  #### D8
    geno.w <- dcast(geno[pond=="D8"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

    geno.w.ag <- geno.w[,.N,chr]
    geno.w <- geno.w[chr%in%geno.w.ag[N>1000]$chr]

    geno.w.comb <- geno.w[,list(AD_REF.LOW = (1-effPA_D8Male1)*effRD_D8Male1 + (1-effPA_D8Male2)*effRD_D8Male2,
                                AD_ALT.LOW = (effPA_D8Male1)*effRD_D8Male1 + (effPA_D8Male2)*effRD_D8Male2,
                                AD_REF.HIGH = (1-effPA_D8PE1)*effRD_D8PE1 + (1-effPA_D8PE2)*effRD_D8PE2,
                                AD_ALT.HIGH = (effPA_D8PE1)*effRD_D8PE1 + (effPA_D8PE2)*effRD_D8PE2),
                            list(CHROM=chr, POS=pos)]

    geno.w.comb[,deltaSNP:=(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - (AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
    geno.w.comb[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
    geno.w.comb[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
    geno.w.comb[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]

    geno.w.comb <- geno.w.comb[REF_FRQ>0.05 & REF_FRQ<0.95]
