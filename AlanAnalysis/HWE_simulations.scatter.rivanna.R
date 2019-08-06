#BiocManager::install("SeqArray")
#install.packages("data.table")
#install.packages("foreach")
#BiocManager::install("SeqVarTools")


  ### libraries
    library(SeqArray)
    library(data.table)
    library(foreach)
    library(doMC)
    registerDoMC(20)
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

  ### generate expected distribution of HWE
      ### expected frequencies
        hwe.stat[,efAA:=afreq^2]
        hwe.stat[,efAa:=2*afreq*(1-afreq)]
        hwe.stat[,efaa:=(1-afreq)^2]

        hwe.stat[,n:=nAA + nAa + naa]

      ### simulation function
        setkey(hwe.stat, n)
        u <- unique(hwe.stat$n)

        simGeno <- function(n.sim=1, sigdig=2, mac=4, args) {
          ##  sigdig <- 2 ### for rounding of genotype frequencies
          ##  mac <- 4 ### minor allele count

          sim.o <- foreach(s=1:n.sim)%do%{

            ### generate random call
             hwe.stat.r <- foreach(i=u)%dopar%{

               print(paste(s, which(i==u), length(u), sep=" / "))

               probs <- as.matrix(hwe.stat[J(i),c("efAA", "efAa", "efaa"), with=F])
               tmp.out <- t(apply(Hmisc::rMultinom(probs, i), 1, function(x) as.numeric(table(factor(x, levels=c("efAA", "efAa", "efaa"))))))
               tmp.out <- as.data.table(tmp.out)
               setnames(tmp.out, names(tmp.out), c("rAA", "rAa", "raa"))

               cbind(hwe.stat[J(i)], tmp.out)
             }

             hwe.stat.r <- rbindlist(hwe.stat.r)
             hwe.stat.r[,frAA:=rAA/(n)]
             hwe.stat.r[,frAa:=rAa/(n)]
             hwe.stat.r[,fraa:=raa/(n)]

           ### calculate normalized genotype density values


             hwe.ag.obs <- hwe.stat.r[afreq>mac/n & afreq<(n-mac)/n,
                                     list(n.obs=length(p), f.hat.obs=mean(f), n=mean(n), nAA=mean(nAA), nAa=mean(nAa), vAa=var(nAa)),
                                     list(fAA=round(fAA, sigdig), fAa=round(fAa, sigdig), faa=round(faa, sigdig),
                                          pond, year)]

             hwe.ag.sim <- hwe.stat.r[afreq>mac/n & afreq<(n-mac)/n,
                                     list(n.sim=length(p), f.hat.sim=mean(f)),
                                     list(fAA=round(frAA, sigdig), fAa=round(frAa, sigdig), faa=round(fraa, sigdig),
                                         pond, year)]

            setkey(hwe.ag.obs, fAA, fAa, faa, pond, year)
            setkey(hwe.ag.sim, fAA, fAa, faa, pond, year)

            hwe.ag.m <- merge(hwe.ag.obs, hwe.ag.sim)
            hwe.ag.m[,en:=log2(n.obs/n.sim)]
            hwe.ag.m[,diff:=n.obs - n.sim]

            hwe.ag.m[,sim:=s]
            hwe.ag.m[,sim.arg:=args]

            hwe.ag.m
          }
          sim.o <- rbindlist(sim.o)
          sim.o
        }

      ### simulate
        hwe.ag.m <- simGeno(n.sim=1, args=args.vec)
        hwe.ag.m[,py:=paste(pond, year, sep=".")]
        save(hwe.ag.m, file=paste("/scratch/aob2x/daphnia_hwe_sims/hwe_sim", args.vec, ".Rdata", sep=""))
