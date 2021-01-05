## ijob -c1 -p standard -A berglandlab
#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

### libraries
  library(data.table)
  library(SeqArray)

### open GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")

### load dpfiltsnps
  load("/project/berglandlab/Karen/MappingDec2019/dpfiltsnps_20200224.Rdata")

### obtusa, pulicaria, and other bad clones
  badClones <- c("Spring_2017_DBunk_340" , "March20_2018_DBunk_39" , "Fall_2016_D10_54" , "March20_2018_D8_19" , "March15_2019_DBunk_MomPE1" ,
                "March15_2019_DBunk_MomPE20" , "April_2017_Dbarb_11" ,"March20_2018_DBunk_26" , "March20_2018_DBunk_37" ,
                "March20_2018_DBunk_42" , "March20_2018_DBunk_10" ,"March20_2018_DBunk_18" , "March20_2018_DBunk_21" ,
                "March20_2018_DBunk_22" , "March20_2018_DBunk_23" ,"March20_2018_DBunk_38" , "March20_2018_DBunk_40" ,
                "March20_2018_DBunk_41" , "March20_2018_DBunk_43" ,"2018_Pulicaria_Pond21_22" , "2018_Pulicaria_Pond22_21" ,
                "2018_Pulicaria_Pond22_53" , "2018_Pulicaria_Pond22_62" ,"2018_Pulicaria_Pond22_72")


### get good samples
  samps <- seqGetData(genofile, "sample.id")
  goodClones <- samps[!samps%in%badClones]

### get allele frequencies
  seqSetFilter(genofile, sample.id=goodClones, variant.id=dpfiltsnps$variant.ids)

  snps.dt <- data.table(variant.ids=seqGetData(genofile, "variant.id"),
                        af=seqAlleleFreq(genofile, .progress=T))

### filter down to ones that are polymorphic in pulex
  pulexPoly <- snps.dt[af!=1 & af!=0]

### merge dpfiltsnps
  setkey(dpfiltsnps, variant.ids)
  setkey(pulexPoly, variant.ids)

  dpfiltsnps <- merge(dpfiltsnps, pulexPoly)
