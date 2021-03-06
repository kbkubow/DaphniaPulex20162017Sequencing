Let's try doing IBS and superclone assignment on the filtered, LD pruned dataset (405 individuals, 155,898 SNPs)
```
#Load libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(foreach)
        library(SeqArray)
        library(doMC)
        registerDoMC(20)
        library(ggplot2)
        library(viridis)
        library(tidyverse)
        
# Load genotype file and snp and sample filter files

        genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
        
        load("finalsetsnpset01_20190423.Rdata")
        
        load("secondsampstokeep_20190422.Rdata")

# Run IBS

        ibs <- snpgdsIBS(genofile, snp.id=finalsetsnpset01, sample.id=secondsampstokeep, num.thread=20, maf=0.01,
          missing.rate=0.15, autosome.only = FALSE)
        
        ### a bit of re-formating of the ibs matrix
                ibs.mat <- ibs$ibs
                rownames(ibs.mat) <- ibs$sample.id
                colnames(ibs.mat) <- ibs$sample.id

  ### quick plot
        ggplot(melt(ibs.mat), aes(Var1,Var2, fill=value)) + geom_raster() + scale_fill_viridis()

        ### make the IBs matrix long form
                ibs.long <- as.data.table(melt(ibs.mat))
                setnames(ibs.long, names(ibs.long), c("cloneA", "cloneB", "distance"))
                ibs.long <- na.omit(ibs.long)
                
        ### Identify cutoff by graphing
        
                ggplot(data=ibs.long, aes(x=distance)) + geom_histogram(binwidth=0.005)               
 ```
 ![InitialIBSDistonFiltSNPsInd](InitialIBSDistonFiltSNPsInd.tiff)

 The presence of multiple D. obtusa is compacting all the pulex at one end of the graph. Let's zoom in on the pulex region.
 
 ```
                 ggplot(data=ibs.long[distance > 0.75], aes(x=distance)) + geom_histogram(binwidth=0.005)
```
 ![InitialIBSDistonFiltSNPsIndZoomA](InitialIBSDistonFiltSNPsIndZoomA.tiff)

 This looks better. There is a long tail between the peak at 1 (individuals within a superclone) and the next peak (comparisons between superclones).
 Let's zoom in a bit more to that region.
 
  ![InitialIBSDistonFiltSNPsIndZoomB](InitialIBSDistonFiltSNPsIndZoomB.tiff)

 Based on this, I am making a decision to make the threshold for superclone inclusion at 0.9775. We can discuss this.

```
  ### Identify superclones
                ### give temporary labels to super clones based on identity of clone B
                        superclone <- ibs.long[distance>.9775]
                        superclone[,SC.sub := as.numeric(as.factor(cloneB))]

                ### collapse nested superclones
                superclone.o <- foreach(sc.i=unique(superclone$SC.sub), .combine="c")%do%{
                        paste(sort(unique(superclone[cloneB%in%superclone[SC.sub==sc.i]$cloneA]$cloneA)), collapse=";")
                }

                superclone.o.uniq <- unique(superclone.o)

                sc.dt <- foreach(i=superclone.o.uniq, .combine="rbind")%do%{
                        data.table(clone=strsplit(i, ";")[[1]],
                                                superClone.size=length(strsplit(i, ";")[[1]]),
                                                superClone.index=which(i==superclone.o.uniq))
                }
                sc.dt[,superClone.sizeRank := as.numeric(as.factor(rank(-superClone.size, ties="average")))]

                sc.dt <- as.data.table(sc.dt %>% 
                        mutate(SCnum = group_indices_(sc.dt, .dots=c("superClone.sizeRank", "superClone.index"))))

    ### do plot to make sure our head is screwed on correctly
      plot(superClone.sizeRank ~ superClone.size, sc.dt)

    ### label superclones with letters. What do you do when you have more than 26 superclones?
        sc.dt[,SC:=LETTERS[SCnum]]
        sc.dt$SC <- ifelse(sc.dt$SCnum==27, "AA", ifelse(sc.dt$SCnum==28, "AB", ifelse(sc.dt$SCnum==29, "AC",
            ifelse(sc.dt$SCnum==30, "AD", ifelse(sc.dt$SCnum==31, "AE", sc.dt$SC)))))

 
    ### rename singletone individuals to "OO" to follow Karen's convention
        sc.dt[superClone.size==1, SC:="OO"]
        
        write.table(sc.dt, file="Superclones20161718_20190423", sep="\t", row.names=FALSE, quote=FALSE)
```
Should also have a superclone file for all individuals, including low coverage, except for two individuals that basically failed.
```
#Load libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(foreach)
        library(SeqArray)
        library(doMC)
        registerDoMC(20)
        library(ggplot2)
        library(viridis)
        library(tidyverse)
        
# Load genotype file and snp and sample filter files

        genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")
        
        load("finalsetsnpset01_20190423.Rdata")
        
        sample.ids <- seqGetData(genofile, "sample.id")
        sampleidsdt <- as.data.table(sample.ids)
        sampleidsnosuperlow <- sampleidsdt$sample.ids[sample.ids!="March20_2018_DBunk_39" & sample.ids!="March20_2018_D8_19" & sample.ids!="Fall_2016_D10_54"]

# Run IBS

        ibs <- snpgdsIBS(genofile, snp.id=finalsetsnpset01, sample.id=sampleidsnosuperlow, num.thread=20, maf=0.01,
          missing.rate=0.15, autosome.only = FALSE)
        
        ### a bit of re-formating of the ibs matrix
                ibs.mat <- ibs$ibs
                rownames(ibs.mat) <- ibs$sample.id
                colnames(ibs.mat) <- ibs$sample.id

  ### quick plot
        ggplot(melt(ibs.mat), aes(Var1,Var2, fill=value)) + geom_raster() + scale_fill_viridis()

        ### make the IBs matrix long form
                ibs.long <- as.data.table(melt(ibs.mat))
                setnames(ibs.long, names(ibs.long), c("cloneA", "cloneB", "distance"))
                ibs.long <- na.omit(ibs.long)
                
        ### Identify cutoff by graphing
        
                ggplot(data=ibs.long, aes(x=distance)) + geom_histogram(binwidth=0.005)               
                ggplot(data=ibs.long[distance > 0.75], aes(x=distance)) + geom_histogram(binwidth=0.005)

  ### Identify superclones
                ### give temporary labels to super clones based on identity of clone B
                        superclone <- ibs.long[distance>.9775]
                        superclone[,SC.sub := as.numeric(as.factor(cloneB))]

                ### collapse nested superclones
                superclone.o <- foreach(sc.i=unique(superclone$SC.sub), .combine="c")%do%{
                        paste(sort(unique(superclone[cloneB%in%superclone[SC.sub==sc.i]$cloneA]$cloneA)), collapse=";")
                }

                superclone.o.uniq <- unique(superclone.o)

                sc.dt <- foreach(i=superclone.o.uniq, .combine="rbind")%do%{
                        data.table(clone=strsplit(i, ";")[[1]],
                                                superClone.size=length(strsplit(i, ";")[[1]]),
                                                superClone.index=which(i==superclone.o.uniq))
                }
                sc.dt[,superClone.sizeRank := as.numeric(as.factor(rank(-superClone.size, ties="average")))]

                sc.dt <- as.data.table(sc.dt %>% 
                        mutate(SCnum = group_indices_(sc.dt, .dots=c("superClone.sizeRank", "superClone.index"))))

    ### do plot to make sure our head is screwed on correctly
      plot(superClone.sizeRank ~ superClone.size, sc.dt)

    ### label superclones with letters. What do you do when you have more than 26 superclones?
        sc.dt[,SC:=LETTERS[SCnum]]
        sc.dt$SC <- ifelse(sc.dt$SCnum==27, "AA", ifelse(sc.dt$SCnum==28, "AB", ifelse(sc.dt$SCnum==29, "AC",
            ifelse(sc.dt$SCnum==30, "AD", ifelse(sc.dt$SCnum==31, "AE", ifelse(sc.dt$SCnum==32, "AF", 
            ifelse(sc.dt$SCnum==33, "AG", ifelse(sc.dt$SCnum==34, "AH", ifelse(sc.dt$SCnum==35, "AI", sc.dt$SC)))))))))

 
    ### rename singletone individuals to "OO" to follow Karen's convention
        sc.dt[superClone.size==1, SC:="OO"]
        
        write.table(sc.dt, file="Superclones20161718withlowcoverageind_20190424", sep="\t", row.names=FALSE, quote=FALSE)
