
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

pat.datll <- rbind(datll[variable.1=="p1"][!variable%in%c("A", "B")],
                    datll[variable%in%c("A", "B")])


mat.datll <- rbind(datll[variable.1=="p2"][!variable%in%c("A", "B")],
                    datll[variable%in%c("A", "B")])












ggplot(data=mat.datll[chr=="Scaffold_9199_HRSCAF_10755"],
        aes(x=as.factor(haplo), y=id, fill=as.factor(value))) +
geom_tile() +
facet_wrap(as.factor(variable=="A" | variable=="B")~chr, scales="free", ncol=1) +
coord_flip()




ggplot(data=pat.datll[chr=="Scaffold_2217_HRSCAF_2652"],
        aes(x=as.factor(haplo), y=id, fill=as.factor(value))) +
geom_tile() +
facet_wrap(as.factor(variable=="A" | variable=="B")~chr, scales="free", ncol=1) +
coord_flip()



ggplot(data=datll,
        aes(x=as.factor(haplo), y=id, fill=as.factor(value))) +
geom_tile() +
facet_grid(variable.1~chr, scales="free") +
coord_flip()
