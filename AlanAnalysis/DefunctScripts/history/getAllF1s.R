#module load intel/18.0 intelmpi/18.0 R/3.6.3; R


library(data.table)

dat <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/D8ParentOffspringRelationships.csv")
setnames(dat, "FocalIndividual", "cloneid")

write.csv(dat[ParentASC=="A" & ParentBSC=="C"][,"cloneid",with=F],
          file="/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.allF1s.delim",
          quote=F, row.names=F)
