#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(SeqArray)

### load posterior probabilities
  setwd("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm")
  loadDat <- function(fn) {
    print(fn)
    #fn <- fns[1]
    pp <- fread(fn, skip="magicReconstruct-Summary,genoprob,Conditonal genotype probability", header=T, fill=T)

    setnames(pp, names(pp), as.character(pp[1,]))
    pp <- pp[grepl("genotype", marker)]

    ppl <- melt(pp[-(1:3)], id.vars="marker")
    ppl[value!="1",value2:=as.numeric(paste(value, "0", sep=""))]
    ppl[value=="1", value2:=1]

    ppl[,genotype:=last(tstrsplit(marker, "_"))]
    ppl[,clone:=gsub("_genotype[0-9]{1,}", "", marker)]

    ppl.ag <- ppl[,list(n=sum(value2>.95)), genotype]

    setkey(ppl, genotype)
    setkey(ppl.ag, genotype)

    ppl <- ppl[J(ppl.ag[n>10]$genotype)]

    ### gneotype code
    gc <- fread(fn, skip="Genotype,Code,founder", nrows=10)
    setnames(gc, "Genotype", "genotype")

    ppl <- merge(ppl, gc)
    codes <- data.table(geno=LETTERS[1:4], founder=names(table(ppl$founder))[c(1,3,2,4)])

    ppl <- merge(ppl, codes, by="founder")
    ppl[,chr:=gsub(".all.out.post.csv", "", fn)]
    return(ppl)
  }
  fns <- list.files("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm", pattern=".out.post.csv")
  ppl <- foreach(x=fns)%do%loadDat(x)
  ppl <- rbindlist(ppl)

  ppl.ag <- ppl[,list(founder=geno[which.max(value2)]), list(variable, clone, chr)]


### make snp.dt
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       key="chr")
  ppl.ag <- merge(ppl.ag, snp.dt, by.x="variable", by.y="id")
  write.csv(ppl.ag, file="/project/berglandlab/alan/AxC_F1.csv", quote=F, row.names=F)

### plot
  plotOut <- ggplot() +
  geom_tile(data=ppl, aes(x=variable, y=founder, fill=value2)) +
  geom_line(data=ppl.ag, aes(x=variable, y=founder, group=clone), color="red") +
  facet_grid(chr.x~clone) +
  scale_fill_viridis() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
