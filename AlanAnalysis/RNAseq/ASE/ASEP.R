library(ASEP)
library(data.table)
library(patchwork)
load("~/ase_geno_phase.star.Rdata")

ase.geno.phase[,id:=paste(chr, pos, sep="_")]
ase.simple <- ase.geno.phase[ref_dosage==1][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")][,c("genes", "samp", "id", "xCount", "totalCount"), with=F]
setnames(ase.simple, c("genes", "samp", "id", "xCount", "totalCount"), c("gene", "id", "snp", "ref", "total"))
ase.simple[,snp:=as.numeric(as.factor(snp))]


ase.gene.table <- ASE_detection(dat_all = ase.simple[gene=="Daphnia00787"], phased=TRUE, varList=NULL, adaptive=TRUE, n_resample=10^3, parallel=FALSE, save_out=FALSE)
dat_M0_phased = as.data.table(phasing(ase.simple[gene=="Daphnia00787"], phased=FALSE, n_condition="one"))

plot_ASE(ase.simple[gene%in%paste("Daphnia0078", c(5:9), sep="")], phased=T)

ggplot(data=ase.simple[gene=="Daphnia00787"], aes(x=id, y=ref/total)) + geom_boxplot()

a1 <- ggplot(data=ase.geno.phase[genes=="Daphnia00787"], aes(x=as.numeric(as.factor(id.x)), y=superclone, fill=as.factor(allele.x))) + geom_tile()
a2 <- ggplot(data=ase.geno.phase[genes=="Daphnia00787"], aes(x=as.numeric(as.factor(id.x)), y=superclone, fill=as.factor(allele.y))) + geom_tile()


ggplot(data=ase.simple[gene%in%paste("Daphnia0078", c(7), sep="")],
              aes(x=id, y=ref/total, fill=clone)) + geom_boxplot() +
              facet_grid(~gene)





load("~/ase_geno_phase.star.Rdata")

ase.simple <- ase.geno.phase[ref_dosage==1][allele.x!=allele.y][superclone=="C"][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")]
ggplot() +
geom_hline( yintercept=mean(ase.simple$xCount/ase.simple$totalCount, na.rm=T)) +
geom_boxplot(data=ase.simple[genes%in%paste("Daphnia0078", c(5:9), sep="")],
              aes(x=clone, y=xCount/totalCount, fill=samp, group=interaction(samp, ref_dosage))) +
facet_grid(~genes, scales="free_x") + theme_bw()

                            ggplot(data=ase.geno.phase[ref_dosage==1][genes%in%paste("Daphnia0078", c(7), sep="")], 
                                          aes(x=id, y=xCount/totalCount, fill=superclone, group=clone)) + geom_boxplot() +
                                          facet_grid(~genes, scales="free_x")





a1 / a2
