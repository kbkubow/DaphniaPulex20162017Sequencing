dt <- data.table(sim.generation = c(1:10000))
dt[,n:=cos((sim.generation - 1) / 100) * 500 + 1000]

plot(n~sim.
generation, dt)


library(ggplot2)
library(data.table)
library(cowplot)


dt <- fread("~/slim/Untitled.txt")
dt[,fracA := (nZW / n)]
dt[,fracB := (nZZ / n)]

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
