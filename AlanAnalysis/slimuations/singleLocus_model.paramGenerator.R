N = c(10000)  #population size
CR = seq(from=.01, to=.99, by=.1) #cloning rate
AFR = seq(from=.9, to=.99, by=.05) #A female rate
BFR = seq(from=.5, to=.99, by=.05) #B female rate


dt <- as.data.table(expand.grid(N, CR, AFR, BFR))
setnames(dt, names(dt), c("N", "CR", "AFR", "BFR"))
dt <- cbind(data.table(id=c(1:dim(dt)[1])), dt)

write.csv(dt, quote=F, row.names=F,
         file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/slimuations/singleLocus_model.paramList")
