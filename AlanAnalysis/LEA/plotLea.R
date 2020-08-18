


# scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/Q_samps.Rdata ~/Q_samps.Rdata
# scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/dap.snmf.Rdata ~/dap.snmf.Rdata

# R

### libraries
  library(LEA)
  library(foreach)
  library(data.table)
  library(ggplot2)
  library(pophelper)
  library(cowplot); theme_set(theme_cowplot())

### load datda
  load("~/Q_samps.Rdata")
  load("~/dap.snmf.Rdata")


### which K?
  plot(dap.snmf)


### extract Q matrix
  processQ <- function(k, samp, orderby) {
    #k<-3; samp=samps$clone; orderby="set"
  	qi<-as.data.table(q.list[[k]])
  	qi <- as.data.table(qi)
  	qi[,gr:=apply(qi, 1, which.max)]
  	qi[,samp:=samp]
  	qil <- melt(qi, measure.vars=paste("V", 1:k, sep=""))
    setnames(qil, "samp", "clone")
    qil <- merge(qil, samps, by="clone")
    qil[set=="pulex_Dcat", set:="pulex_DCat"]



    if(orderby=="gr") {
    	qil.ag <- foreach(gr.i=unique(qil$gr), .combine="rbind")%do%{
    		qil.sub <- qil[gr==gr.i][variable==paste("V", gr.i, sep="")]
    		qil.sub[,ord:=rank(value, ties.method="first")]
    		qil.sub
    	}
    } else if(orderby=="set") {
      qil.ag <- foreach(set.i=unique(qil$set), .combine="rbind")%do%{
        #set.i<-"pulex_D8"
        qil.temp <- qil[set==set.i]
        qil.temp.ag <- qil.temp[,list(gr=unique(gr), max=max(value)), list(clone)]
        keycol <-c("gr", "max")
        setorderv(qil.temp.ag, keycol)
        qil.temp.ag[,ord:=c(1:dim(qil.temp.ag)[1])]
        qil.temp.ag
      }
    }

  	setkey(qil, clone)
  	setkey(qil.ag, clone)

  	qil <- merge(qil, qil.ag)
  	qil[,ordf:=factor(ord, levels=sort(unique(ord)))]
    qil[,set:=factor(set,
                    levels=c("obtusa_DBunk",
                            "pulicaria_Pond22",
                            "pulex_W6",
                            "pulex_W1",
                            "pulex_D10",
                            "pulex_D8",
                            "pulex_DBunk",
                            "pulex_DCat",
                            "pulex_DOil",
                            "pulex_Dramp"))]
  	return(qil)
  }

  plotQ <- function(k.i) {

    ggplot(k.i, aes(x=ord, y=value, fill=as.factor(variable))) +
    geom_bar(position="stack", stat="identity", width = 1) +
    facet_grid(~set, scale="free_x", space="free") +
    geom_text(data=k[Species=="pulicaria"], aes(x=ord, y=1, label="*"), size=2) +
    geom_text(data=k[Species=="obtusa"], aes(x=ord, y=1, label="$"), size=2) +
    theme(strip.text.x = element_text(angle=90),
          panel.spacing = unit(.2, "lines"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_x_continuous(expand = c(0,0)) +
    xlab("")

  }


plot.q1 <- plotQ(processQ(k=1, samp=samps$clone, orderby="set"))
plot.q2 <- plotQ(processQ(k=2, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q3 <- plotQ(processQ(k=3, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q4 <- plotQ(processQ(k=4, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q5 <- plotQ(processQ(k=5, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q6 <- plotQ(processQ(k=6, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q7 <- plotQ(processQ(k=7, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())

plot.q1 / plot.q2 / plot.q3 / plot.q4 / plot.q5 / plot.q6 / plot.q7



hist(k[variable.x=="V6"]$value.x)
k[variable.x=="V6"][value.x>.2]


hist(log10(k[variable.x=="V3"]$value.x), breaks=1000)
hist(log10(k[variable.x=="V6"]$value.x), breaks=1000)
