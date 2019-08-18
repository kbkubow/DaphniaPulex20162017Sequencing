### libraries
  library(data.table)
  library(cowplot)
  library(ggplot2)

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


        ≈
        theme_bw


  /Users/alanbergland/Documents/GitHub/
