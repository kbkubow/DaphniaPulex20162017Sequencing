
### ported from DaphniaPulex20162017Sequencing/Filter201620172018VCFinR/IBSSuperCloneonFilteredSnpSet.md


### libraries
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

### Load genotype file and snp and sample filter files

  genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

  load("/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/finalsnpstousewSimoids_20190430.Rdata")

  load("/mnt/spicy_3/Karen/201620172018FinalMapping/ForAlan/secondsampstokeepwSimo_20190430.Rdata")

# Run IBS

    ibs <- snpgdsIBS(genofile, snp.id=finalsnpstousewSimoids, sample.id=secondsampstokeep, num.thread=20, maf=0.01,
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

    ### The presence of multiple D. obtusa is compacting all the pulex at one end of the graph. Let's zoom in on the pulex region.

  ggplot(data=ibs.long[distance > 0.75], aes(x=distance)) + geom_histogram(binwidth=0.005)

  ### This looks better. There is a long tail between the peak at 1 (individuals within a superclone) and the next peak (comparisons between superclones).
  ### Let's zoom in a bit more to that region.


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
        sc.labeler <- function(scnum) {
          mod <- scnum - 26
          if(mod<=0) return(LETTERS[scnum])
          if(mod>0) return(paste("A", LETTERS[scnum-26], sep=""))
        }
        sc.labeler <- Vectorize(sc.labeler)

        sc.dt[superClone.size>1, SC := sc.labeler(SCnum)]

        sc.dt[superClone.size==1, SC:="OO"]

########## plot

  ### sort IBS matrix to make it look pretty (the code looks cluky as hell, but it works)

    ### first, kick out obtusa
      sc <- fread("/mnt/spicy_3/AlanDaphnia/vcf/Superclones20161718updated_20190802.csv")
      sc.dt <- sc.dt[clone%in%sc[Species=="Pulex"]$clone]

    ### remake ibs.long
        ibs.mat.fig <- ibs.mat
        ibs.long <- as.data.table(melt(ibs.mat.fig))
        setnames(ibs.long, names(ibs.long), c("cloneA", "cloneB", "distance"))
        ibs.long <- na.omit(ibs.long)

    ### first, need to tack in SC identities
        setnames(ibs.long, "cloneA", "clone")
        setkey(ibs.long, "clone")
        setkey(sc.dt, "clone")
        ibs.long <- merge(ibs.long, sc.dt)
        setnames(ibs.long, "clone", "cloneA")
        setnames(ibs.long, "SC", "SC.A")

        setnames(ibs.long, "cloneB", "clone")
        setkey(ibs.long, "clone")
        setkey(sc.dt, "clone")
        ibs.long <- merge(ibs.long, sc.dt)
        setnames(ibs.long, "clone", "cloneB")
        setnames(ibs.long, "SC", "SC.B")


	### re-rank SCs based on size in pond

        ibs.long <- ibs.long[,c("cloneA", "cloneB", "SC.A", "SC.B", "distance")]
		ibs.long[,pondA := tstrsplit(cloneA, "_")[[3]]]
   	 	ibs.long[,pondB := tstrsplit(cloneB, "_")[[3]]]


		ibs.long.ag <- ibs.long[,list(pond.n = length(distance)), list(pondA, SC.A) ]
		ibs.long.ag[,pond.sc.rank := ibs.long.ag[,list(pond.sc.rank = rank(-pond.n, ties="random")), list(pondA)]$pond.sc.rank]
		ibs.long.ag[,pond.sc.rank := letters[pond.sc.rank]]

		### be lazy and write a loop
			ibs.long.2 <- foreach(i=1:dim(ibs.long.ag)[1], .combine="rbind")%do%{
					temp <- ibs.long[pondA==ibs.long.ag$pondA[i] & SC.A==ibs.long.ag$SC.A[i]]
					temp[,sc.a:=ibs.long.ag$pond.sc.rank[i]]
					temp
			}

			ibs.long.3 <- foreach(i=1:dim(ibs.long.ag)[1], .combine="rbind")%do%{
					temp <- ibs.long.2[pondB==ibs.long.ag$pondA[i] & SC.B==ibs.long.ag$SC.A[i]]
					temp[,sc.b:=ibs.long.ag$pond.sc.rank[i]]
					temp
			}

			ibs.long <- ibs.long.3

	### next, generate [s]uper[c]lone[i]ds for individual 'A' and 'B'

        ibs.long[,scid.a := paste(SC.A, sprintf("%03d", as.numeric(as.factor(cloneA))), sep=".")]
        ibs.long[,scid.b := paste(SC.B, sprintf("%03d", as.numeric(as.factor(cloneB))), sep=".")]

    ### group on pond

        ibs.long[,pondA := factor(pondA, levels=c("D10", "D8", "DBunk", "DCat", "DLily", "DOil", "Dcat", "Dramp", "W1", "W6"))]
        ibs.long[,pondB := factor(pondB, levels=c("D10", "D8", "DBunk", "DCat", "DLily", "DOil", "Dcat", "Dramp", "W1", "W6"))]

        ibs.long[,scid.a := paste(LETTERS[as.numeric(pondA)], scid.a, sep=".")]
        ibs.long[,scid.b := paste(LETTERS[as.numeric(pondB)], scid.b, sep=".")]

        ibs.long[,scid.a := as.factor(scid.a)]
        ibs.long[,scid.b := as.factor(scid.b)]

    ### tack in buffer cloneIds for graphical purposes
        ibs.long <- rbind(ibs.long,
                          data.table(scid.a=c(paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "DCat", "W1")]$pondA)], min(ibs.long$SC.A), sep=".")), c("000"), sep="."),
                                              paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "Dramp", "W6")]$pondA)], max(ibs.long$SC.A), sep=".")), c("999"), sep=".")),
                                     scid.b=c(paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "DCat", "W1")]$pondB)], min(ibs.long$SC.B), sep=".")), c("000"), sep="."),
                                              paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "Dramp", "W6")]$pondB)], max(ibs.long$SC.B), sep=".")), c("999"), sep="."))),
                          fill=T)
         ibs.long[,scid.a := factor(scid.a, levels=sort(unique(as.character(scid.a))))]
         ibs.long[,scid.b := factor(scid.b, levels=sort(unique(as.character(scid.b))))]

    ### make lower triangle poofy-de-poof
        ibs.long[,dist.noTri := distance]
        ibs.long[as.numeric(scid.a)>as.numeric(scid.b), dist.noTri:=NA]

    ### make pond bounding boxes
        ibs.long.ag <- data.table(scid.a.min=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "DCat", "W1")]$pondA)], min(ibs.long$SC.A, na.rm=T), sep=".")), c("000"), sep="."),
                                  scid.a.max=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "Dramp", "W6")]$pondA)], max(ibs.long$SC.A, na.rm=T), sep=".")), c("999"), sep="."),
                                  scid.b.min=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "DCat", "W1")]$pondB)], min(ibs.long$SC.B, na.rm=T), sep=".")), c("000"), sep="."),
                                  scid.b.max=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "Dramp", "W6")]$pondB)], max(ibs.long$SC.B, na.rm=T), sep=".")), c("999"), sep="."))

        ibs.long.ag[,scid.a.min := as.numeric(factor(scid.a.min, levels=sort(unique(as.character(ibs.long$scid.a)))))]
        ibs.long.ag[,scid.a.max := as.numeric(factor(scid.a.max, levels=sort(unique(as.character(ibs.long$scid.a)))))]
        ibs.long.ag[,scid.b.min := as.numeric(factor(scid.b.min, levels=sort(unique(as.character(ibs.long$scid.b)))))]
        ibs.long.ag[,scid.b.max := as.numeric(factor(scid.b.max, levels=sort(unique(as.character(ibs.long$scid.b)))))]

    ### save
      save(ibs.long, ibs.long.ag, file="~/ibs.Rdata")

    ### plot it
        h.just <- .25
        v.just <- .25
        l.size <- 1.5

        ibs.long[SC.A=="A" & SC.B=="A", x.lab:="A"]
        ibs.long[SC.A=="B" & SC.B=="B", y.lab:="B"]

      pdf("~/ibs.pdf", h=8, w=9)
       ggplot(data=ibs.long, aes(x=scid.a, y=scid.b, fill=distance)) +
        geom_raster() +
        scale_fill_viridis(option="D")
      dev.off()


        +
        geom_rect(xmin=ibs.long.ag$scid.a.min[1]-2*h.just, xmax=ibs.long.ag$scid.a.max[1]+h.just,
                  ymin=ibs.long.ag$scid.b.min[1]-2*v.just, ymax=ibs.long.ag$scid.b.max[1]+v.just,
                  fill=NA, color="red", size=l.size) +
        geom_rect(xmin=ibs.long.ag$scid.a.min[2]-2*h.just, xmax=ibs.long.ag$scid.a.max[2]+h.just,
                  ymin=ibs.long.ag$scid.b.min[2]-2*v.just, ymax=ibs.long.ag$scid.b.max[2]+v.just,
                  fill=NA, color="red", size=l.size) +
        geom_rect(xmin=ibs.long.ag$scid.a.min[3]-2*h.just, xmax=ibs.long.ag$scid.a.max[3]+h.just,
                  ymin=ibs.long.ag$scid.b.min[3]-2*v.just, ymax=ibs.long.ag$scid.b.max[3]+v.just,
                  fill=NA, color="red", size=l.size) +
        geom_rect(xmin=ibs.long.ag$scid.a.min[4]-2*h.just, xmax=ibs.long.ag$scid.a.max[4]+h.just,
                  ymin=ibs.long.ag$scid.b.min[4]-2*v.just, ymax=ibs.long.ag$scid.b.max[4]+v.just,
                  fill=NA, color="red", size=l.size) +
        geom_rect(xmin=ibs.long.ag$scid.a.min[5]-2*h.just, xmax=ibs.long.ag$scid.a.max[5]+h.just,
                  ymin=ibs.long.ag$scid.b.min[5]-2*v.just, ymax=ibs.long.ag$scid.b.max[5]+v.just,
                  fill=NA, color="red", size=l.size)








































        #write.table(sc.dt, file="Superclones20161718_20190423", sep="\t", row.names=FALSE, quote=FALSE)
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
