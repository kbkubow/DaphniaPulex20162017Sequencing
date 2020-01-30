### libraries
 library(gdsfmt)
 library(SNPRelate)
 library(data.table)
 library(ggplot2)
 library(foreach)
 library(lattice)
 library(tidyr)
 library(SeqArray)


## open files
  genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapDec19PulexOnly_filtsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

## Load SNPs, filtered but not LD pruned.
  load("snpsvarpulexpresentinhalf_20200114.Rdata")
## Load SNPs, fixedAB
  load("ABfixed_20200115.Rdata")
  ABfixedids <- ABfixed$variant.ids
  seqSetFilter(genofile, variant.id=ABfixedids)


  seqSetFilter(genofile, variant.id=snpsvarpulexpresentinhalf)

### Load individuals, let's try first 2 As and 2 Bs and 2018 F1s

  load("F1sandAB.Rdata")
  F1sandABids <- as.character(F1sandAB$clone)

  seqSetFilter(genofile, sample.id=F1sandABids)

  D82016cloneid <- c("Spring_2016_D8_8.12", "Spring_2016_D8_8.3", "Spring_2016_D8_8.18")
  seqSetFilter(genofile, sample.id=D82016cloneid)

   snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
         chr = seqGetData(genofile, "chromosome"),
         pos = seqGetData(genofile, "position"),
         dp = seqGetData(genofile, "annotation/info/DP"))

   snps9200 <- snps[chr=="Scaffold_1863_HRSCAF_2081"]

   snps9200input <- snps9200$variant.ids

     seqSetFilter(genofile, variant.id=snps9200input)


    ### get residual heterozygosity

     het <- t(seqGetData(genofile, "$dosage"))
     het <- as.data.table(het)

     colnames(het) <- c(seqGetData(genofile, "sample.id"))
     het$variant.ids <- seqGetData(genofile, "variant.id")

     setkey(het, variant.ids)
     setkey(snps, variant.ids)

     mhet <- merge(snps, het)

     #mhethet <- mhet[May_2017_D8_515==0 & May_2017_D8_612==1 |
          May_2017_D8_515==0 & May_2017_D8_612==2 |
          May_2017_D8_515==1 & May_2017_D8_612==0 |
          May_2017_D8_515==1 & May_2017_D8_612==2 |
          May_2017_D8_515==2 & May_2017_D8_612==0 |
          May_2017_D8_515==2 & May_2017_D8_612==1 |
          May_2017_D8_515==1 & May_2017_D8_612==1 ]

     #mhet <- mhethet

     #mhet[, chr:=as.integer(chr)]


     #sample.ids <- seqGetData(genofile, "sample.id")
     sample.ids <- F1sandABids

     normdosagetable <- foreach(i=sample.ids, .combine="rbind")%do%{

     data.table(chr=mhet$chr, pos=mhet$pos, clone=c(i), normdosage=ifelse(mhet[[i]]==2 & mhet$Spring_2016_D8_8.12==0, 0,
       ifelse(mhet[[i]]==0 & mhet$Spring_2016_D8_8.12==0, 2, mhet[[i]])))

     }

     numclones <- data.table(clone=c(sample.ids), numclone=c(1:25))

     setkey(normdosagetable, clone)
     setkey(numclones, clone)

     normdosagetablenum <- merge(normdosagetable, numclones)

     numpos <- data.table(pos=c(mhet$pos), numpos=c(1:49428))

     setkey(normdosagetablenum, pos)
     setkey(numpos, pos)

     normdosagetablenumpos <- merge(normdosagetablenum, numpos)

     ggplot(data=normdosagetablenumpos, aes(x=numpos, y=numclone, color=as.factor(normdosage))) + geom_point()
     ggplot(data=normdosagetablenumpos[numpos < 5000], aes(x=numpos, y=numclone, color=as.factor(normdosage))) + geom_point()
     ggplot(data=normdosagetablenumpos[numpos > 5000 & numpos < 10000], aes(x=numpos, y=numclone, color=as.factor(normdosage))) + geom_point()
     ggplot(data=normdosagetablenumpos[numpos > 10000 & numpos < 15000], aes(x=numpos, y=numclone, color=as.factor(normdosage))) + geom_point()
     ggplot(data=normdosagetablenumpos[numpos < 20000], aes(x=numpos, y=numclone, color=as.factor(normdosage))) + geom_point()


   ### Adding in superclones

   superclone01 <- fread("superclonemaf01.csv")

   colnames(superclone01) <- c("clone", "supercloneloose", "superclonecons")

   setkey(normdosagetablenumpos, clone)
   setkey(superclone01, clone)

   msc <- merge(normdosagetablenumpos, superclone01)

   msc$numclone <- factor(msc$numclone, levels = msc$numclone[order(msc$supercloneloose)])

   msc[, normdosage:=as.factor(normdosage)]


   ggplot(data=msc[numpos < 2000], aes(x=numpos, y=numclone, color=normdosage)) + geom_point() + theme(text = element_text(size=7))
   ggplot(data=msc[numpos > 250 & numpos < 1250], aes(x=numpos, y=numclone, color=normdosage)) + geom_point() + theme(text = element_text(size=7))


11500 13300

   setkey(msc, chr, pos)
   setkey(ABFKfixsub, chr, pos)

   mscwinfo <- merge(msc, ABFKfixsub)

   mscwinfo$BFKreads <- ifelse(mscwinfo$homrefA=="0", mscwinfo$refreads, mscwinfo$altreads)
   mscwinfo$Areads <- ifelse(mscwinfo$homrefA=="22", mscwinfo$refreads, mscwinfo$altreads)
   mscwinfo$propB <- mscwinfo$BFKreads/(mscwinfo$BFKreads+mscwinfo$Areads)
