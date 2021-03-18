scp aob2x@rivanna.hpc.virginia.edu:~/cdlo_250K.boot.trees.Rdata ~/.

### libraries
  library(ggplot2)
  library(data.table)
  library(ggtree)
  library(ape)
  library(foreach)

### data
  load("~/cdlo_250K.boot.trees.Rdata")

  cd.tree.name <- foreach(i=1:length(cdl.tree))%do%{
    message(i)
    tree.i<-cdl.tree[i]
    w <- tstrsplit(tree.i[[1]]$tip.label[1], "fa_")[[2]]
    w
  }
  head(cd.tree.name)
  cd.tree.name <- unlist(cd.tree.name)

  names(cdl.tree) <- cd.tree.name

  cd.tree.name[grepl("Scaffold_9197_HRSCAF_10753", cd.tree.name)]

njo <- cdl.tree$"Scaffold_9197_HRSCAF_10753_2452444-2702443"
njo <- root(njo, "pulicaria.1.fa_Scaffold_9197_HRSCAF_10753_2452444-2702443")


d <- data.table(label=njo$tip.label)
d[,label:=gsub("Dcat", "DCat", label)]
d[,pond:=tstrsplit(label, "_")[[3]]]


d[,A:=ifelse(grepl("April_2017_D8_213", label), "A", "")]
d[,C:=ifelse(grepl("April_2017_D8_151", label), "C", "")]

d[grepl("March20_2018_Dcat_2", label)]

d[A!="" | C!="" | pond=="DCat" | pond=="Dcat", clean:=label]
d[grepl("March20_2018_Dcat_2", label)]

d[is.na(clean), clean:=""]
d <- as.data.frame(d)

ggtree(njo) %<+% d +
geom_tiplab(aes(label=A)) +
geom_tiplab(aes(label=C)) +
geom_tiplab(aes(label=clean, color=pond), offset=0, angle=0, align=F)  +
geom_tiplab(aes(label=pond, color=pond), offset=0, angle=0, align=T)  +
theme(legend.position = "none")+ ggtitle("QTL10")
