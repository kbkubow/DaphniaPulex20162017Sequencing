
library("data.table")
library("ggplot2")

dat <- fread("/mnt/sammas_storage/bergland-lab/alan/testTrio.consensus.header.phase.csv")
setnames(dat, c("#CHROM", "POS"), c("chr", "pos"))
setkey(dat, chr, pos)


dat <- dat
dat[,id:=c(1:dim(dat)[1])]
datl <- melt(dat, id.var=c("chr", "pos", "REF", "ALT", "id"))
datl[,p1:=substr(value, 0, 1)]
datl[,p2:=substr(value, 3, 3)]

datll <- melt(datl, id.var=c("chr", "pos", "REF", "ALT", "id", "variable"), measure.var=c("p1", "p2"))

datll[,haplo:=paste(variable, variable.1, sep="_")]


ggplot(data=datll[chr=="Scaffold_7757_HRSCAF_8726"][pos>=(8660157-5000) & pos<=(8660157+5000)],
        aes(x=as.factor(haplo), y=id, fill=as.factor(value))) +
geom_tile() +
facet_wrap(~chr, scales="free") +
coord_flip()



ggplot(data=datll[chr=="Scaffold_7757_HRSCAF_8726"],
        aes(x=as.factor(haplo), y=id, fill=as.factor(value))) +
geom_tile() +
facet_wrap(variable.1~chr, scales="free") +
coord_flip()


ggplot(data=datll,
        aes(x=as.factor(haplo), y=id, fill=as.factor(value))) +
geom_tile() +
facet_grid(variable.1~chr, scales="free") +
coord_flip()
