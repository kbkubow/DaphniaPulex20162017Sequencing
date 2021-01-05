### libraries
  library(data.table)
  library(cowplot)
  library(ggplot2)

### mesocosm asex/sex output data
  dat <- fread("DaphniaPulex20162017Sequencing/AlanAnalysis/DataFiles/SampleCountsMeso2019Week6.csv")


  dat.ag <- dat[,list(TotalInd=mean(TotalInd, na.rm=T),
                      Embryos=mean(EstimatedEmbryosTot, na.rm=T),
                      Male=mean(PropTotMale, na.rm=T)),
                 list(Week, SC, Clone)]
  dat.ag[is.nan(Embryos), Embryos:=0]
  dat.ag[is.nan(Male), Male:=0]

  popsize.plot <- ggplot(data=dat.ag[SC!="AB"][Week<=6], aes(x=Week, y=TotalInd, group=Clone, color=SC)) +
  geom_line(size=1.5) + ylab("Individuals / Liter") +
  theme(legend.position="none")

  embryos.plot <- ggplot(data=dat.ag[SC!="AB"][Week<=6], aes(x=Week, y=Embryos, group=Clone, color=SC)) +
  geom_line(size=1.5) + ylab("Num. sexual emryos") +
  theme(legend.position="none")

  male.plot <- ggplot(data=dat.ag[SC!="AB"][Week<=6], aes(x=Week, y=Male, group=Clone, color=SC)) +
  geom_line(size=1.5) + ylab("Proportion Male") +
  theme(legend.position="none")


  ggsave(plot_grid(popsize.plot, embryos.plot, male.plot, nrow=1), file="~/popGrowth.pdf", height=4, width=9)

#### Hatchling data
  ### load data
    load("DaphniaPulex20162017Sequencing/AlanAnalysis/DataFiles/hatchsummary20190814.Rdata")

  ### Survival
    summary(glm(cbind(PropSurv*SurvN, (1-PropSurv)*SurvN) ~ Vernalized + Clone,
                data=hatchsummary[Clone%in%c("AxB", "D8222", "D8515")],
                family=binomial()))

    surv.plot <- ggplot(data=hatchsummary[Clone%in%c("AxB", "D8222", "D8515")],
          aes(x=c("Not Vern", "Vern")[Vernalized+1], y=PropSurv, color=as.factor(Vernalized))) +
          facet_grid(~Clone) +
          geom_errorbar(aes(x=c("Not Vern", "Vern")[Vernalized+1],
                            ymin=PropSurv - 1.96*sqrt((PropSurv*(1-PropSurv)/SurvN)),
                            ymax=PropSurv + 1.96*sqrt((PropSurv*(1-PropSurv)/SurvN))),
                        width=.0) +
          geom_point(size=4) +
          theme(legend.position="none") +
          xlab("") + ylab("Survival Rate") + panel_border(colour="black", size=.5)

  ### Fertility
    summary(glm(cbind(PropFertile*FertN, (1-PropFertile)*FertN) ~ Vernalized + Clone,
                data=hatchsummary[Clone%in%c("AxB", "D8222", "D8515")],
                family=binomial()))

    fert.plot <- ggplot(data=hatchsummary[Clone%in%c("AxB", "D8222", "D8515")],
          aes(x=c("Not Vern", "Vern")[Vernalized+1], y=PropFertile, color=as.factor(Vernalized))) +
          facet_grid(~Clone) +
          geom_errorbar(aes(x=c("Not Vern", "Vern")[Vernalized+1],
                            ymin=PropFertile - 1.96*sqrt((PropFertile*(1-PropFertile)/FertN)),
                            ymax=PropFertile + 1.96*sqrt((PropFertile*(1-PropFertile)/FertN))),
                        width=.0) +
          geom_point(size=4) +
          theme(legend.position="none") +
          xlab("") + ylab("Fertility Rate") + panel_border(colour="black", size=.5)

  ### Male induction
    summary(glm(cbind(PropMaleObs*MaleN, (1-PropMaleObs)*MaleN) ~ Clone,
                data=hatchsummary[Clone%in%c("AxB", "D8222", "D8515")],
                family=binomial()))

    male.plot <- ggplot(data=hatchsummary[Clone%in%c("AxB", "D8222", "D8515")],
          aes(x=c("Not Vern", "Vern")[Vernalized+1], y=PropMaleObs, color=as.factor(Vernalized))) +
          facet_grid(~Clone) +
          geom_errorbar(aes(x=c("Not Vern", "Vern")[Vernalized+1],
                            ymin=PropMaleObs - 1.96*sqrt((PropMaleObs*(1-PropMaleObs)/MaleN)),
                            ymax=PropMaleObs + 1.96*sqrt((PropMaleObs*(1-PropMaleObs)/MaleN))),
                        width=.0) +
          geom_point(size=4) +
          theme(legend.position="none") +
          xlab("") + ylab("Male Rate") + panel_border(colour="black", size=.5)


    #### joint plot
      plot_grid(surv.plot, fert.plot, male.plot, ncol=1)
