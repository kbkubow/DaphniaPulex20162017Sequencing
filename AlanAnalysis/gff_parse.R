# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(foreach)
  library(data.table)

### gff
  gff <- fread("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff")
  #gff[,ID:=gsub("ID=", "", first(tstrsplit(tstrsplit(V9, ";|:")[[1]], "-")))]
  gff[,ID:=gsub("ID=", "", tstrsplit(V9, ";|:")[[1]])]

  gff[,ss:=paste(V1, V3, V4, V5, ID, sep="_")]
  setkey(gff, ss)

  gff2 <- gff[!duplicated(gff)]
  dim(gff)
  dim(gff2)

  collapse_fun <- function(start, stop) {
    #start <- c(1,5); stop<-c(4, 10)
    length(unique(unlist(apply(cbind(start, stop), 1, function(x) x[1]:x[2]))))
  }


  gff.ag <- gff2[V3=="CDS",list(cds=sum(V5-V4 + 1), cds2=collapse_fun(start=V4, stop=V5)), list(chr=V1, id=ID)]
  gff.ag[,aa:=cds2/3]
  gff.ag
  summary(gff.ag)
  gff.ag[aa!=round(aa)]

  gff.ag[,gene:=first(tstrsplit(id, "-"))]


### use this
  gff.ag.ag <- gff.ag[,list(max_cds=max(cds)), list(gene)]
