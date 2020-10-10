# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  #library(qtl)
  library(data.table)
  #library(ggplot2)
  #library(patchwork)
  library(lme4)

############################
#### Prep the input data ###
############################

########################
### [P]henotype data ###
########################
  ### Load SC data
    sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020//Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  ## Male rates
    ### load raw male phenotype data
      males <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/MaleCensus.csv")
      #males <- fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rqtl.convertData.R")

    ### some renaming
     sconeliter <- sc[OneLiterPheno==1]
     sconeliter$SCB <- ifelse(sconeliter$LabGenerated==0 & sconeliter$SC!="OO", sconeliter$SC, sconeliter$clone)
     sconeliter$Clone <- sconeliter$clone

    ### tidy
     setkey(males, Clone)
     setkey(sconeliter, Clone)
     mmales <- merge(males, sconeliter)
     mmales$propmale <- mmales$Males/mmales$NewTotal
     mmales$propmalenoneo <- mmales$Males/(mmales$NewTotal-mmales$Neos)

     male <- mmales[,c("Clone", "Replicate", "Males", "NewTotal", "propmale", "SCB"), with=F]

     ### tack in whether it has been genotyped yet
       #f1s.use <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.onlyPheno.delim")
       #f1s.use <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.allF1s.delim")

       f1s <- sc[AxCF1Hybrid==1][OneLiterPheno==1]$clone
       f1s.use <- data.table(cloneid=f1s)
       setnames(f1s.use, "cloneid", "Clone")
       f1s.use[,geno:=T]

       setkey(mmales, "Clone")
       setkey(f1s.use, "Clone")
       mmales <- merge(mmales, f1s.use, all.x=T, all.y=T)

  ## Epphipia fill rates and production
    ### load raw epphiphial fill data
     epp <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/EphippiaFinal.csv")
     sconeliter <- sc[OneLiterPheno==1]
     sconeliter$SCB <- ifelse(sconeliter$LabGenerated==0 & sconeliter$SC!="OO", sconeliter$SC, sconeliter$clone)
     sconeliter$Clone <- sconeliter$clone

     setkey(epp, Clone)
     setkey(sconeliter, Clone)
     mepp <- merge(epp, sconeliter)
     epp <- mepp

   # Add new variables
     epp$TotalEppB <- epp$TotalEpp*2
     epp$fill <- epp$TotalEmbCorr/epp$TotalEppB
     epp$fullid <- paste(epp$Clone,epp$Rep,sep="_")
     epp$Rep <- as.factor(epp$Rep)
     epp$Clone <- as.factor(epp$Clone)
     epp$Dayssince1stDate <- as.factor(epp$Dayssince1stDate)
     epp$SCB <- as.factor(epp$SCB)

     setkey(epp, "Clone")
     setkey(f1s.use, "Clone")
     epp <- merge(epp, f1s.use, all.x=T, all.y=T)


### summaries and merging
  ### male
    mmales.ag <- mmales[,list(nClones=length(unique(Clone)),
                              anySeq=any(geno),
                              clone=Clone[geno==T][!is.na(Clone)][1],
                              propmale=sum(Males)/sum(NewTotal),
                              propmalenoneo=sum(Males)/sum(NewTotal-Neos),
                              nMales=mean(Males),
                              nDaps=mean(NewTotal),
                              nNeo=mean(Neos)),
                         list(SCB,
                              population, AxCF1Hybrid)][!is.na(SCB)]


    epp.ag <- epp[,list(nClones=length(unique(Clone)),
                              anySeq=any(geno),
                              fill=mean(fill, na.rm=T),
                              fill.se=sd(fill, na.rm=T)/sqrt(sum(TotalEppB)),
                              epp=mean(TotalEppB, na.rm=T)),
                         list(SCB, population, AxCF1Hybrid)][!is.na(SCB)]



### averages and BLUPs
  ### averages
      epp.ag[SCB=="A", gr:="A"]
      epp.ag[SCB=="C", gr:="C"]
      epp.ag[(grepl("AxB", SCB) | AxCF1Hybrid==1) & anySeq==T, gr:="AxC"]
      epp.ag[grepl("AxB", SCB) & is.na(anySeq), gr:="CxC"]

      setkey(mmales.ag, SCB, population, anySeq)
      setkey(epp.ag, SCB, population, anySeq)


    ### generate BLUPs
    male.mod <- glmer(propmalenoneo ~ 1 + (1|SCB) + (1|Clone:Replicate),
                      data=mmales,
                      family=binomial(),
                      weights=NewTotal-Neos)

    fill.mod <- glmer(fill~1+(1|Dayssince1stDate)+(1|Clone:Rep)+(1|SCB),
                data=epp,
                family=binomial(),
                weights=TotalEppB)

    epp.mod <- glmer(TotalEppB~1+(1|Dayssince1stDate)+(1|Clone:Rep)+(1|SCB),
                data=epp,
                family=poisson())


    r1 <- data.table(SCB=rownames(ranef(male.mod)$SCB),
                    propmalenoneo.ranef =  ranef(male.mod)$SCB[,1])
    r2 <- data.table(SCB=rownames(ranef(fill.mod)$SCB),
                    fill.ranef =  ranef(fill.mod)$SCB[,1])
    r3 <- data.table(SCB=rownames(ranef(epp.mod)$SCB),
               epp.ranef =  ranef(epp.mod)$SCB[,1])

    r <- merge(r1, r2, by="SCB", all.x=T, all.y=T)
    r <- merge(r, r3, by="SCB", all.x=T, all.y=T)


    m <- merge(mmales.ag, epp.ag, all.x=T, all.y=T)

    mm <- merge(m, r, by="SCB")


### save
  save(mm, r, mmales, mmales.ag, epp, epp.ag, file="~/F1_pheno.Rdata")

### load pheno & geno data

### libraries
  library(data.table)
  library(lme4)
  library(doMC)
  registerDoMC(10)
  library(glmperm)

  ### Load SC data
    sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020//Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  ### load phenotype data
    load(file="~/F1_pheno.Rdata")

  ### load genotype data (from DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.convert_output.vitPath.R)
    pio <- fread(file= "/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all_AxC.vitPath.pio.csv")
    pio[,phase.fold.geno:=phase.geno]
    pio[phase.fold.geno%in%c(12, 21), phase.fold.geno:=12]

    ((table(pio$imputedGeno, pio$phase.fold.geno)))

### simple gwas
pio.ag <- pio[,list(.N), list(chr, pos, id)]

ii <- seq(from=1, to=dim(pio.ag)[1], by=1)


out <- foreach(i=ii, .errorhandling="remove")%dopar%{
  print(paste(which(i==ii), length(ii), sep=" / "))

  # id.i <- 2142150 - max Gprime peak
  # i <- 5158
  # id.i <- 3065470
  id.i <- pio.ag$id[i]

  ### male model: phase.geno
    tmp.male <- merge(mmales, pio[id==id.i], by="clone")
    tmp.male.ag <- tmp.male[,list(Males=sum(Males), N=sum(NewTotal - Neos), phase.fold.geno=phase.fold.geno[1], ng=length(unique(phase.fold.geno))), list(SCB)]

    t1.male <- glmer(propmalenoneo ~ as.factor(phase.fold.geno) + (1|SCB) + (1|Replicate:SCB), tmp.male, family=binomial(),
                    weights=tmp.male$NewTotal - tmp.male$Neos)
    t0.male <- glmer(propmalenoneo ~ 1 + (1|SCB) + (1|Replicate:SCB), tmp.male, family=binomial(),
                    weights=tmp.male$NewTotal - tmp.male$Neos)

    t1.male <- glm(cbind(Males, N-Males) ~ as.factor(phase.fold.geno), tmp.male.ag, family=binomial())
    t0.male <- glm(cbind(Males, N-Males) ~ 1, tmp.male.ag, family=binomial())

    aov.male <- anova(t1.male, t0.male, test="Chisq")

  ### fill rate model
    tmp.fill <- merge(epp, pio[id==id.i], by="clone")
    tmp.fill.ag <- tmp.fill[,list(emb=sum(fill*TotalEppB, na.rm=T), N=sum(TotalEppB), phase.fold.geno=phase.fold.geno[1], ng=length(unique(phase.fold.geno))), list(SCB)]
    tmp.fill.ag[,perm:=sample(phase.fold.geno)]



    t1.fill <- glmer(fill ~ as.factor(phase.fold.geno) + (1|SCB) + (1|Rep:SCB), tmp.fill, family=binomial(),
                    weights=tmp.fill$TotalEppB)
    t0.fill <- glmer(fill ~ 1 + (1|SCB) + (1|SCB:Rep), tmp.fill, family=binomial(),
                    weights=tmp.fill$TotalEppB)
    anova(t1.fill, t0.fill, test="Chisq")
    
    t1.fill <- glm(cbind(emb, N-emb) ~ as.factor(phase.fold.geno) , tmp.fill.ag, family=quasibinomial())
    t0.fill <- glm(cbind(emb, N-emb) ~ 1 , tmp.fill.ag, family=quasibinomial())
    aov.fill <- anova(t1.fill, t0.fill, test="Chisq")

    tmp.fill.ag[,perm:=sample(phase.fold.geno)]
    t1.fill <- glm(cbind(emb, N-emb) ~ as.factor(perm) , tmp.fill.ag, family=binomial())
    t0.fill <- glm(cbind(emb, N-emb) ~ 1 , tmp.fill.ag, family=binomial())
    anova(t1.fill, t0.fill, test="Chisq")


  ### output
    #data.table(id=id.i,
    #            chisq=c(aov.fill[2,6], aov.male[2,6]),
    #            p=c(aov.fill[2,8], aov.male[2,8]),
    #            term=c("fill", "male"))
    averages <- merge(tmp.male.ag[,list(males.mu=sum(Males)/sum(N), males.n=sum(N), geno.n=length(N)), list(phase.fold.geno)],
                      tmp.fill.ag[,list(fill.mu=sum(emb)/sum(N), fill.n=sum(N), geno.n=length(N)), list(phase.fold.geno)],
                      by="phase.fold.geno")
    averages[,id:=id.i]

    out <- data.table(id=id.i,
                chisq=c(aov.fill[2,4], aov.male[2,4]),
                p=c(aov.fill[2,5], aov.male[2,5]),
                term=c("fill", "male"))
    merge(out, averages, by="id", allow.cartesian=T)

}

out <- rbindlist(out)
out[,pa:=p.adjust(p, "fdr")]
setkey(out, id)
setkey(pio.ag, id)

out <- merge(out, pio.ag)
save(out, file="~/gwasOutput.Rdata")
out[term=="fill"][which.min(chisq)]
out[id==977045]

out[id==3065515]





scp aob2x@rivanna.hpc.virginia.edu:~/gwasOutput.Rdata ~/.

library(data.table)
library(ggplot2)
load("~/gwasOutput.Rdata")
load("~/peaks.Rdata")


id.i <- sapply(c(1:12), function(x) out[chr==peaks$CHROM[x]][which.min(abs(peaks$posMaxGprime[x] - pos))]$id)

id.i <- c(out[chr==peaks$CHROM[11]][which.min(abs(peaks$posMaxGprime[11] - pos))]$id,
          out[chr==peaks$CHROM[11]][which.min(abs(peaks$start[11] - pos))]$id,
          out[chr==peaks$CHROM[11]][which.min(abs(peaks$end[11] - pos))]$id)

out[id==id.i[6]]


ggplot() +
geom_line(data=out, aes(x=id, y=-log10(p))) +
geom_point(data=out[id%in%id.i], aes(x=id, y=-log10(p)), color="red") +
facet_wrap(~term)


eppresid.plot <- ggplot() +
#geom_hline(yintercept=summary(perm)[2]) +
geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
geom_rect(data=peaks, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="grey", alpha=.5) +
geom_line(data=mr.epp.dt, aes(x=pos, y=lod, color=chr), size=1) +
facet_grid(.~chr, scales="free_x") +
theme(strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
ggtitle("Num. Epp")


#setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))
Gprime.plot <- ggplot(data=gprime, aes(x=pos, y=Gprime, color=chr)) +
#geom_rect(data=peaks, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="grey", alpha=.5) +
geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
geom_line(size=.75) +
facet_grid(.~chr, scales="free_x") +
theme(strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
ggtitle("pooledWild")









ggplot(data=out, aes(x=id, y=-log10(p))) + geom_point() + facet_wrap(~term)

### test
load(file="/nv/vol186/bergland-lab/alan/peaks.Rdata")

pio[chr==peaks$CHROM[11]][which.min(abs(peaks$posMaxGprime[11] - pos))]
tmp.male <- merge(mmales, pio[id==3042479], by="clone")


t1.male <- glmer(propmalenoneo ~ as.factor(phase.geno) + (1|SCB), tmp.male, family=binomial(),
                  weights=tmp.male$NewTotal - tmp.male$Neos)
summary(t1.male)


anova(lm(propmalenoneo~as.factor(phase.geno), tmp))
t1.male <- lm(propmalenoneo~as.factor(imputedGeno), tmp)
