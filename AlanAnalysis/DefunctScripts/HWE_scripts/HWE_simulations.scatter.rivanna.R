#BiocManager::install("SeqArray")
#install.packages("data.table")
#install.packages("foreach")
#BiocManager::install("SeqVarTools")


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

  ### load precomputed HWE data
    load(file="/scratch/aob2x/daphnia_hwe_sims/hwe_stat.Rdata")

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
             hwe.stat.r <- foreach(i=u)%do%{

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
        write.csv(hwe.ag.m, quote=F, row.names=F,
                  file=paste("/scratch/aob2x/daphnia_hwe_sims/hweSimOut/hwe_sim", args.vec, ".csv", sep=""))
