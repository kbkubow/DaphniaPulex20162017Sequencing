#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### opebn genofiles
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")

### load hap data
  haps.files <- list.files(path="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/", pattern="*haps")

  haps <- foreach(fi=haps.files, .combine="rbind")%do%{
    #fi <- haps.files[1]
    message(fi)
    tmp <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/", fi, sep="/"))
    tmp[,cloneid:=tstrsplit(V1, " ")[[1]]]
    tmp[,chr:=tstrsplit(V1, " ")[[2]]]
    tmp[,allele1:=tstrsplit(V4, "\\|")[[1]]]
    tmp[,allele2:=tstrsplit(V4, "\\|")[[2]]]

    tmp
  }
  haps[V4=="A_m|C_m", encode:=1]
  haps[V4=="A_p|C_m", encode:=2]
  haps[V4=="A_m|C_p", encode:=3]
  haps[V4=="A_p|C_p", encode:=4]









  seqSetFilter(genofile, unique(haps$chr))

  dt <- data.table(V2=seqGetData(genofile, "variant.id"), startPos=seqGetData(genofile, "position"))
  haps <- merge(haps, dt, by="V2")

  haps[,allele1:=V4]
  haps[,allele2:=V5]

  haps[V4=="A1", allele1:="A"]
  haps[V4=="A2", allele1:="B"]
  haps[V4=="C2", allele1:="D"]
  haps[V4=="C1", allele1:="C"]

  haps[V5=="A1", allele2:="A"]
  haps[V5=="A2", allele2:="B"]
  haps[V5=="C2", allele2:="D"]
  haps[V5=="C1", allele2:="C"]

  haps[,diplo:=apply(cbind(haps$allele1, haps$allele2), 1, function(x) paste(sort(x), collapse=""))]

  prop.table(table(haps$V1, haps$diplo))



  sum(haps[is.na(encode)]$stopPos - haps[is.na(encode)]$startPos)
  sum(haps[!is.na(encode)]$stopPos - haps[!is.na(encode)]$startPos)
