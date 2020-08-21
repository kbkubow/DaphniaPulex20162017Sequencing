


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
  #load("~/Q_samps.Rdata")
  #load("~/dap.snmf.Rdata")
  load("~/snmf_out.Rdata")

### which K?
  ggplot(data=ce.dt, aes(x=k, y=ce, group=run, color=as.factor(run))) + geom_line()
  ce.dt[k<10][which.min(ce)]

### extract Q matrix
  processQ <- function(k, samp, run, orderby) {
    #k<-8; run=8; samp=samps$clone; orderby="set"
  	qi<-as.data.table(q.list[[k]][[run]])
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
        qil.temp.ave <- qil.temp[,list(mu=mean(value)), list(variable)]


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
    geom_text(data=k.i[SC=="A"], aes(x=ord, y=1.1, label="A"), size=5) +
    geom_text(data=k.i[SC=="C"], aes(x=ord, y=1.1, label="C"), size=5) +

    theme(strip.text.x = element_text(angle=90),
          panel.spacing = unit(.2, "lines"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_x_continuous(expand = c(0,0)) +
    xlab("") +
    guides(fill=FALSE)

  }

  plotQ.sp <- function(k.i) {
    k.i[SC=="A",col:="A"]
    k.i[SC=="C",col:="C"]
    setorderv(k.i, "col")

    puli.k <- k.i[Species=="pulicaria"][which.max(value)]$variable
    obt.k <- k.i[Species=="obtusa"][which.max(value)]$variable

    puli.plot <- ggplot(data=k.i[variable==puli.k], aes(x=set, y=log10(value), color=col)) +
    geom_jitter(width = 0.25) +
    #geom_jitter(data=k.i[variable==puli.k][SC=="A"], aes(x=set, y=log10(value)), color="blue", size=1, width = 0.25) +
    #geom_jitter(data=k.i[variable==puli.k][SC=="C"], aes(x=set, y=log10(value)), color="red", size=1, width = 0.25) +
    coord_flip() +
    labs(title=puli.k) +
    theme(legend.position="none")

    obt.plot <- ggplot(data=k.i[variable==obt.k], aes(x=set, y=log10(value), color=col)) +
    geom_jitter(width = 0.25) +
    #geom_jitter(data=k.i[variable==puli.k][SC=="A"], aes(x=set, y=log10(value)), color="blue", size=1, width = 0.25) +
    #geom_jitter(data=k.i[variable==puli.k][SC=="C"], aes(x=set, y=log10(value)), color="red", size=1, width = 0.25) +
    coord_flip() +
    labs(title=obt.k) +
    theme(legend.position="none")


    puli.plot | obt.plot

  }



### summary plot
  ce.dt[k==8][which.min(ce)]

  k.i<-8; run.i<-8
  kr <- processQ(k=k.i, run=run.i, samp=samps$clone, orderby="set")
  q.plot <- plotQ(kr)
  out.plot <- plotQ.sp(kr)
  ce.plot <-   ggplot(data=ce.dt, aes(x=k, y=ce, group=run, color=as.factor(run))) +
               geom_line() +
               geom_point(data=ce.dt[k==k.i][run==run.i], aes(x=k, y=ce), size=3, color="black") +
               theme(legend.position="none")


layout <- "
AAAAAAAA
##BBBBCC
"
  q.plot / ( out.plot | ce.plot) +
  plot_layout(design = layout)


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
