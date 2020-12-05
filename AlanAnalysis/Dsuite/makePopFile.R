
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(SeqArray)

### open GDS file
  #genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")

### clone metadata file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  sc <- sc[Nonindependent==0]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]

### hard filtering of SC
  sc.ag <- sc[LabGenerated==F, list(clone=clone[which.max(medrd)][1]), list(SC.uniq, Species)]
  #sc.ag[,clone:=paste(Species, clone, sep="_")]


sc.ag[,pond:=tstrsplit(clone, "_")[[3]]]
sc.ag[,gr:=pond]
sc.ag[Species!="pulex", gr:=Species]
sc.ag[gr=="obtusa", gr:="Outgroup"]

### write to disk
  write.table(sc.ag[,c("clone", "gr"),with=F], file="/scratch/aob2x/daphnia_hwe_sims/dsuite/SETS.txt", quote=F, row.names=F, col.names=F, sep="\t")




#### write to disk
#  write.table(sc.ag, file="/scratch/aob2x/daphnia_hwe_sims/popPhase/representativeSC.delim", quote=F, row.names=F, sep=",")
#
#
#### generate bcf file lists for merge
#  ### get chrs
#    chrs <- system("cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' | sort | uniq | grep -v 'chr'", intern=T)
#
#    foreach(chr.i=chrs)%do%{
#      #chr.i <- 2
#      message(chr.i)
#      bcfs <- system(paste("ls /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/*", chr.i, ".phase.vcf.gz", sep=""), intern=T)
#
#      bcfs.use <- foreach(clone.i = sc.ag$clone, .combine="c")%do%{
#        #clone.i<-sc.ag$clone[122]
#        grep(paste(clone.i, ".Scaffold", sep=""), bcfs)
#      }
#
#
#      o <- bcfs[bcfs.use]
#      write.table(o, file=paste("/scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/", chr.i, ".list", sep=""), quote=F, row.names=F, col.names=F)
#    }
#
