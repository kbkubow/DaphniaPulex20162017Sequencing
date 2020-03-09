### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)

### load pre-computed file (probably from here: DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/harp_windows.R)
  load("scratch/aob2x/daphnia_hwe_sims/harp_pools/summarizedOut/harpWins.Rdata")

### anyting jump out?
  o[,sex:=gsub("1|2", "", pool)]

  o.ag <- o[,list(freq=c(mean(pA1), mean(pA2), mean(pB1), mean(pB2)),
                  allele=c("A1", "A2", "B1", "B2")),
             list(chr, start, stop, sex)]

  o.ag.ag <- o.ag[,list(delta=(freq[sex=="D8Male"]) - (freq[sex=="D8PE"])),
                   list(chr, start, stop, allele)]

  m <- o.ag.ag[,list(dist=sum(delta^2)), list(chr, start, stop)]

  setkey(m, chr, start, stop)
  setkey(o.ag, chr, start, stop)
  setkey(o.ag.ag, chr, start, stop)

  summary(m)

  save(o.ag, o.ag.ag, m, file="/nv/vol186/bergland-lab/alan/harp_summarized_play.Rdata")


### load
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(foreach)
  library(doMC)
  registerDoMC(4)

  load("/mnt/sammas_storage/bergland-lab/alan/harp_summarized_play.Rdata")
  load(file="/mnt/sammas_storage/bergland-lab/alan/peaks.Rdata")


    o.ag.ag <- o.ag[,list(delta=(freq[sex=="D8Male"]) - (freq[sex=="D8PE"])),
                     list(chr, start, stop, allele)]

   m <- o.ag.ag[,list(dist=sum(delta^2)), list(chr, start, stop)]

   m[,i:=c(1:dim(m)[1])]



   plotPeak <- function(i, buffer=5000) {
    #i<-6; buffer<-5000
     tmp1 <- o.ag[chr==peaks[i]$CHROM][start>=(peaks[i]$start-buffer)][stop<=(peaks[i]$end+buffer)]
     tmp2 <- o.ag.ag[chr==peaks[i]$CHROM][start>=(peaks[i]$start-buffer)][stop<=(peaks[i]$end+buffer)]
     tmp3 <- m[chr==peaks[i]$CHROM][start>=(peaks[i]$start-buffer)][stop<=(peaks[i]$end+buffer)]

     deltaPlot <- ggplot(data=tmp2, aes(x=(start/2+stop/2), y=delta, group=allele, color=allele, linetype=as.factor(grepl("A", allele)))) + geom_line()
     freqPlot <- ggplot(data=tmp1, aes(x=(start/2+stop/2), y=freq, group=interaction(allele, sex), color=sex, linetype=allele)) + geom_line()
     distPlot <- ggplot(data=tmp3, aes(x=(start/2+stop/2), y=dist)) + geom_line() +
                 geom_hline(yintercept=quantile(m$dist, .975, na.rm=T)) + geom_hline(yintercept=quantile(m$dist, .025, na.rm=T))

     plot_grid(freqPlot, deltaPlot, distPlot)
   }

   plotPeak(4)


   ### merge QTL & peaks files
    ### effect plot
       findPeak <- function(i) {
          #i<-1
          setkey(m, chr)
          tmp <- as.data.frame(m)
          tmp <- as.data.table(tmp[m$chr==peaks$chr[i],])
          tmp$posDist <- sqrt((tmp$start - peaks[i]$start)^2 + (tmp$stop - peaks[i]$end)^2)
          tmp <- tmp[which.min(posDist)]

          return(tmp$i)
        }

        peaks[,i:=unlist(sapply(c(1:dim(peaks)[1]), findPeak))]


      ### QTL plot
        setkey(gprime, CHROM, POS)
        gprime.i <- foreach(i=m$i)%dopar%{
          print(paste(i, max(m$i)))
          tmp <- gprime[J(data.table(CHROM=m$chr[i], POS=m$start[i]:m$stop[i])), nomatch=0]
          tmp[,i:=(seq(from=(i-1), to=i, length.out=dim(tmp)[1]) + 0.5) ]
          return(tmp)
        }
        gprime.i <- rbindlist(gprime.i)
        setnames(gprime.i, "CHROM", "chr")


      #### individuals QTL effects
        #tmp.ag <- tmp[,list(delta.mu=mean(delta), delta.sd=sd(delta)), list(allele, sex)]

        ## version 1
          i<-7
          tmp <- o.ag[chr==peaks[i]$chr][start>=peaks[i]$start][stop<=peaks[i]$end]
          tmp.ag <- tmp[,list(freq.mu=mean(freq), freq.sd=sd(freq)), list(allele, sex)]
          ggplot() +
          geom_bar(data=tmp.ag, aes(x=allele, y=freq.mu, group=sex, fill=sex), stat="identity", position="dodge")

        ## version 2
          i<-12
          buffer <- 10000
          tmp <- o.ag.ag[chr==peaks[i]$chr][start>=(peaks[i]$start-buffer)][stop<=(peaks[i]$end+buffer)]
          tmp.ag <- tmp[,list(delta.mu=median(delta), delta.sd=sd(delta)), list(allele)]
          ggplot() +
          geom_bar(data=tmp.ag, aes(x=allele, y=delta.mu, group=allele, fill=allele), stat="identity", position="dodge") +
          geom_errorbar(data=tmp.ag, aes(x=allele, ymin=delta.mu-delta.sd, ymax=delta.mu+delta.sd, group=allele), position=position_dodge(1), width=.1) +
          coord_flip() + guides(fill = guide_legend(reverse = TRUE))




        ### combined

        qtl.plot <- ggplot(data=gprime.i, aes(x=i, y=negLog10Pval, color=chr)) +
        geom_vline(data=peaks, aes(xintercept=i), color="black", linetype="dashed") +
        geom_line(size=1) +
        theme(legend.position = "none")

        effect.plot <- ggplot(data=m, aes(x=i, y=-100*dist, color=chr)) +
        geom_vline(data=peaks, aes(xintercept=i), color="black", linetype="dashed") +
        geom_line(size=.75) +
        theme(legend.position = "none") +
        scale_x_continuous(position = "top")


        plot_grid(qtl.plot, effect.plot, nrow=2)





  ggplot(data=m, aes(x=i, y=dist, color=chr)) + geom_line()

  o.ag.ag[J(m[grepl("7757", chr)][which.max(dist)])]

  o.ag[J(m[grepl("7757", chr)][which.max(dist)])]
