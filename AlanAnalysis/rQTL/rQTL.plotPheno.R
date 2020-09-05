

### a few plot
  # scp aob2x@rivanna.hpc.virginia.edu:~/F1_pheno.Rdata ~/.
  # R
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(cowplot)
  theme_set(theme_cowplot())

  load("~/F1_pheno.Rdata")

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

plot(propmalenoneo~propmalenoneo.ranef, mm)
plot(fill~fill.ranef, mm)
plot(epp~epp.ranef, mm)


  ### weighted means
    fill.rate <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=fill, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05)) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Ephippial Fill Rate")

    num.epp <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=epp, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05)) +
    theme(legend.position = "none")+
    xlab("") +
    ylab("Number of Ephippia")

    propmale <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=propmalenoneo, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05)) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Proportion Male")

  ### ranef
    fill.rate.r <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=fill.ranef, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05)) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Ephippial Fill Rate, BLUP")

    num.epp.r <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=epp.ranef, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05)) +
    theme(legend.position = "none")+
    xlab("") +
    ylab("Number of Ephippia, BLUP")

    propmale.r <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=propmalenoneo.ranef, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05)) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Proportion Male, BLUP")


    (fill.rate | num.epp | propmale) /
    (fill.rate.r| num.epp.r | propmale.r)



### princomp
  pr <- prcomp(na.omit(mm[,c("fill", "epp", "propmalenoneo"), with=F]), center=T, scale=T)

  ggplot(data=mm[!is.na(gr)], aes(y=fill, x=propmalenoneo, color=as.factor(gr))) +
  geom_point()
