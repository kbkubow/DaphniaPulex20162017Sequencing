#module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(SeqArray)

### make SNP table
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile))
  setkey(snp.dt, id)
