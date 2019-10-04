###module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0


### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(foreach)

### load output data
  fl <- system("ls /scratch/aob2x/daphnia_hwe_sims/slim_output/slim_*", intern=T)

  dt <- foreach(i=fl, .combine="rbind")%do%fread(i)

  save(dt, file="~/dt.Rdata")

  load("dt.Rdata")

### some summaries
  dt.ag <- dt[gen>75, list(CF=median(frac_Clonal_Female), CM=median(frac_Clonal_Male),
                          SF=median(frac_Sexual_Female), SM=median(frac_Sexual_Male)),
                    list(N=param_N, CR=param_CR, AFR=param_AFR, BFR=param_BFR)]

  ggplot(data=dt.ag, aes(x=CR, y=CF+CM)) + geom_point() + geom_abline(intercept = 0, slope = 1)


  dt.ag <- dt[gen>50, list(fracA=median(fracA)),
                    list(N=param_N, CR=param_CR, AFR=param_AFR, BFR=param_BFR)]

  ggplot(data=dt.ag, aes(x=CR, y=fracA, color=as.factor(AFR - BFR))) + geom_point() + geom_abline(intercept = 0, slope = 1)







  ggplot(data=dt, aes(x=gen, y=frac_Clonal_Female, group=param_N, color=as.factor(param_N))) + geom_line() + geom_point()


  dt[,fracA := (nZW / n)]
  dt[,fracB := (nZZ / n)]

### summaires
  dt.ag <- dt[gen>10,list(scFrac=c(median(fracA), median(fracB)), sc=c("A", "B"),
                    sexFrac=c(median(frac_Clonal_Female), median(frac_Clonal_Male))), list(param_N)]





dt <- fread("~/slim/Untitled.txt")



dt.long <- melt(dt, id.vars="gen", measure=c("frac_Clonal_Female", "frac_Clonal_Male", "frac_Sexual_Female", "frac_Sexual_Male", "fracA", "fracB"))
dt.long[,sex:=tstrsplit(variable, "_")[[3]]]
dt.long[,mode:=tstrsplit(variable, "_")[[2]]]

popSize.plot <- ggplot(data=dt, aes(x=gen, y=n, color=pop)) +
                  geom_vline(xintercept=c(20, 40, 60, 80), size=.05) + geom_line()

sex.plot <-  ggplot(data=dt.long[!variable%in%c("fracA", "fracB")],
                    aes(x=gen, y=value, group=variable, color=sex, linetype=mode)) +
              geom_vline(xintercept=c(20, 40, 60, 80), size=.05) +
              geom_line() + ylim(0,1)

AB.plot <- ggplot(data=dt.long[variable%in%c("fracA", "fracB")],
                  aes(x=gen, y=value, group=variable, color=variable)) +
            geom_vline(xintercept=c(20, 40, 60, 80), size=.05) +
            geom_line() + ylim(0,1)


plot_grid(popSize.plot, sex.plot, AB.plot, nrow=3, align = "v")
