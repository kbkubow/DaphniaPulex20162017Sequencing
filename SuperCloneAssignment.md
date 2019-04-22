#############################
### Superclone assignment ###
### (originally from an   ###
### email from Karen on   ###
### Feb 8, 2019 "Re:assigning
### super clonse script"  ###
#############################

### Karen edited this on April 22nd and put the script here on GitHub

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

        ### try the phased file
                seqVCF2GDS("/mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.Scaffold_$
                                        "/mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFil$


        ### open gds file
                #genofile <- seqOpen("/mnt/spicy_3/Karen/Sp2017/NewMapping//totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss_lowDpmis$
       # genofile <- seqOpen("/mnt/spicy_3/Karen/Sp2017/NewMapping//totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.seq.gds")
                genofile <- seqOpen("/mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.$


        ### load SNP filtering file
                #load(file="/mnt/spicy_3/AlanDaphnia/outputData/snps.Rdata")

    ### load sample filtering file
        #load(file="/mnt/spicy_3/AlanDaphnia/outputData/sampleids.Rdata")

    ### make IBS matrix
        ### note, there will be some "monomorphic" SNPs in here because we are using a more strictly filtered set of genotype calls (get detai$

                #snp.id=snps[use==T][filter.Dpu.fixed==F]$variant.ids,
                #sample.id=sampleids[goodclone==T]$clone,

                ibs <- snpgdsIBS(genofile,

                                                num.thread=20,
                        maf=0.05,
                        missing.rate=0.25,
                                                autosome.only = FALSE)

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

        ### Identify superclones
                ### give temporary labels to super clones based on identity of clone B
                        superclone <- ibs.long[distance>.975]
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

    ### There was a slight error in assigning letters to the superclones, Karen is adding in a line of code and editing one other line. THIS IS NOT WORKING EITHER!!!!
                sc.dt <- as.data.table(sc.dt %>% 
                        mutate(SCnum = group_indices_(sc.dt, .dots=c("superClone.sizeRank", "superClone.index"))))

    ### do plot to make sure our head is screwed on correctly
      plot(superClone.sizeRank ~ superClone.size, sc.dt)

    ### label superclones with letters
      #sc.dt[,SC := LETTERS[superClone.sizeRank]]
        sc.dt[,SC:=LETTERS[SCnum]]

    ### rename singletone individuals to "O" to follow Karen's convention
      sc.dt[superClone.size==1, SC:="O"]

  ### sort IBS matrix to make it look pretty (the code looks cluky as hell, but it works)
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

        ibs.long <- ibs.long[,c("cloneA", "cloneB", "SC.A", "SC.B", "distance"), with=F]
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
                          data.table(scid.a=c(paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "DCat", "W1")]$p,
                                              paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "Dramp", "W6")]$
                                     scid.b=c(paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "DCat", "W1")]$p$
                                              paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "Dramp", "W6")]$$
                          fill=T)
         ibs.long[,scid.a := factor(scid.a, levels=sort(unique(as.character(scid.a))))]
         ibs.long[,scid.b := factor(scid.b, levels=sort(unique(as.character(scid.b))))]

    ### make lower triangle poofy-de-poof
        ibs.long[,dist.noTri := distance]
        ibs.long[as.numeric(scid.a)>as.numeric(scid.b), dist.noTri:=NA]

    ### make pond bounding boxes
        ibs.long.ag <- data.table(scid.a.min=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "DCat", "W1")]$po$
                                  scid.a.max=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "Dramp", "W6")]$p$
                                  scid.b.min=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "DCat", "W1")]$po$
                                  scid.b.max=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "Dramp", "W6")]$p$

        ibs.long.ag[,scid.a.min := as.numeric(factor(scid.a.min, levels=sort(unique(as.character(ibs.long$scid.a)))))]
        ibs.long.ag[,scid.a.max := as.numeric(factor(scid.a.max, levels=sort(unique(as.character(ibs.long$scid.a)))))]
        ibs.long.ag[,scid.b.min := as.numeric(factor(scid.b.min, levels=sort(unique(as.character(ibs.long$scid.b)))))]
        ibs.long.ag[,scid.b.max := as.numeric(factor(scid.b.max, levels=sort(unique(as.character(ibs.long$scid.b)))))]


    ### plot it
        h.just <- .25
        v.just <- .25
        l.size <- 1.5
       ggplot(data=ibs.long, aes(scid.a, scid.b, fill=distance)) +
        geom_raster() +
        scale_fill_viridis(option="D") +
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


### playing
        library(tsne)
        ibs.mat.tsne <- tsne(ibs.mat)
        ibs.mat.tsne <- as.data.table(ibs.mat.tsne)
        ibs.mat.tsne[,clone:=row.names(ibs.mat)]
        setkey(ibs.mat.tsne, clone)
        ibs.mat.tsne <- merge(ibs.mat.tsne, sc.dt)
        ibs.mat.tsne[,pond := tstrsplit(clone, "_")[[3]]]

        ggplot(data=ibs.mat.tsne[pond=="D8"], aes(V1, V2, color=SC)) + geom_point() + facet_wrap(~pond)

### save
    save(sc.dt, file="/mnt/spicy_3/AlanDaphnia/outputData/superclone.Rdata")
