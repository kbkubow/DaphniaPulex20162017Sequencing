
### librarioes
  library(data.table)
  library(foreach)

### get job ID
  args <- commandArgs(trailingOnly = TRUE)

  nJobs <- as.numeric(args[1])
  jobId <- as.numeric(args[2])

  print(nJobs)
  prin(jobId)
  
#### load data
#  load("/project/berglandlab/Karen/SingleMomsMales20182019/PooledMomsMales/totrdfilt.Rdata")
#  totrdfilt[,effRD:=floor(RRD)]
#  totrdfilt[,effPA:=round(propalt*effRD)/effRD]
#
#
#### make job array table
#  snps.dt <- totrdfilt[,list(pos=unique(pos)), list(chr)]
#  snps.dt[,id:=c(1:dim(snps.dt)[1])]
#  groups <-  split(snps.dt$id, ceiling(seq_along(snps.dt$id)/(length(snps.dt$id)/nJobs)))
#
#### get tmp data
#  setkey(snps.dt, id)
#  snps.dt.tmp <- snps.dt[J(groups[[jobId]])]
#
#  setkey(totrdfilt, chr, pos)
#  setkey(snps.dt.tmp, chr, pos)
#  #totrdfilt.tmp <- totrdfilt[J(snps.dt.tmp)]
#
#### iterate through SNPs
#  o <- foreach(i=1:dim(snps.dt.tmp)[1], .combine="rbind", .errorhandling="remove")%do%{
#    #i<-1
#    tmp <- totrdfilt[J(snps.dt.tmp[i])]
#
#    t1 <- glm(effPA~sex + pond, data=tmp, weights=effRD, family=binomial())
#    t1.D8 <- glm(effPA~sex , data=tmp[pond=="D8"], weights=effRD, family=binomial())
#    t1.DBunk <- glm(effPA~sex , data=tmp[pond=="DBunk"], weights=effRD, family=binomial())
#    t1.DCat <- glm(effPA~sex , data=tmp[pond=="DCat"], weights=effRD, family=binomial())
#
#    data.table(chr=tmp$chr[1], pos=tmp$pos[1],
#               p.value=c(summary(t1)$coef[2,4],
#                         summary(t1.D8)$coef[2,4],
#                         summary(t1.DBunk)$coef[2,4],
#                         summary(t1.DCat)$coef[2,4]),
#                mod=c("all", "D8", "DBunk", "DCat"))
#
#  }
#
#### save
#  write.csv(o, file=paste("/scratch/aob2x/daphnia_hwe_sims/mf/mf_glm_", jobId, ".csv"))
#
#
#
### move to another place
#  fl <- system("ls XYZ", intern=T)
##  o <- foreach(fl.i=fl)%do%{
###    fread(fl.i)
# # }
#  o <- rbindlist(o)
