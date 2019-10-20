###module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0


### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(foreach)
  library(viridis)

### load output data
  fl <- paste("/scratch/aob2x/daphnia_hwe_sims/slim_output/", system("ls /scratch/aob2x/daphnia_hwe_sims/slim_output/", intern=T), sep="")

  dt <- foreach(i=fl, .combine="rbind", .errorhandling="remove")%do%{
    print(which(fl==i))
    #if(which(fl==i)%%25==0) paste(print(which(fl==i), length(fl), sep=" / "))
    fread(i)
  }

  save(dt, file="~/dt.Rdata")

  load("~/dt.Rdata")
  setnames(dt, names(dt), c("pop", "gen", "n",
        "frac_Clonal_Female", "frac_Clonal_Male",
        "frac_Sexual_Female", "frac_Sexual_Male",
        "fZ",
        "nZZ",
        "nZW",
        "nWW", "neut",
        "simID", "N", "CR", "AMR", "BMR"))
    dt[,fZW := (nZW / n)]
    dt[,fZZ := (nZZ / n)]
    #dt[,fZ:= (nZZ+.5*nZW)/n]

### some summaries

  dt.ag <- dt[gen>200, list(fZW=mean(fZW), fZZ=mean(fZZ), fZ=mean(fZ),
                          CF=mean(frac_Clonal_Female), CM=mean(frac_Clonal_Male),
                          SF=mean(frac_Sexual_Female), SM=mean(frac_Sexual_Male),
                          neutFixation=any(neut==0 | neut==1),
                          neutDelta=neut[gen==500] - neut[gen==201],
                          fZdelta=fZ[gen==500] - fZ[gen==201]),
                    list(N=N, CR=CR, maleRate_ZW=AMR, maleRate_ZZ=BMR)]


  eq.freq.plot <- ggplot(data=dt.ag, aes(x=maleRate_ZW, y=maleRate_ZZ, fill=(CM+SM))) +
  geom_raster() + scale_fill_viridis(limits=c(0,1)) + facet_wrap(~CR) + geom_abline(intercept=0, slope=1)  +
  geom_point(data= dt.ag[(1-maleRate_ZZ)==.86 & (1-maleRate_ZW)==.96], aes(x=maleRate_ZW, y=maleRate_ZZ), color="red") +
  geom_point(data= dt.ag[fZZ>=(2*N-1)/(2*N) | fZZ<=(1/(2*N))], aes(x=maleRate_ZW, y=maleRate_ZZ), color="blue", size=.5)

  eq.freq.plot



  ggplot(data=dt.ag[CR<.05 | CR>.99 | CR==.51], aes(x=maleRate_ZW, y=maleRate_ZZ, fill=(fZ))) +
  geom_raster() + scale_fill_viridis(limits=c(0,1)) + facet_wrap(~CR) + geom_abline(intercept=0, slope=1)


   +
  geom_point(data= dt.ag[(1-maleRate_ZZ)==.86 & (1-maleRate_ZW)==.96], aes(x=maleRate_ZW, y=maleRate_ZZ), color="red") +
  geom_point(data= dt.ag[fracA>=(2*N-1)/(2*N) | fracA<=(1/(2*N))], aes(x=maleRate_ZW, y=maleRate_ZZ), color="blue", size=.5)




  summary(lm(fracB~CR*maleRate_ZW*maleRate_ZZ, dt.ag))


  freqB.plot <- ggplot







#### some tests
  library(markovchain)
  library(data.table)
  library(ggplot2)

  dat <- fread("~/slim_out.delim")
  setnames(dat, names(dat), c("pop", "gen", "n",
				"frac_Clonal_Female", "frac_Clonal_Male",
				"frac_Sexual_Female", "frac_Sexual_Male",
				"fZ",
				"nZZ",
				"nZW",
				"nWW", "neut",
				"simID", "N", "CR", "AMR", "BMR"))
  dat[,fracA := (nZW / n)]
  dat[,fracB := (nZZ / n)]
  dat[,fZ:=(nZZ+.5*nZW)/n]
  dat[,seedF:=as.numeric(as.factor(simID))]

  ggplot() +
  geom_line(data=dat, aes(x=gen, y=neut, group=seedF), color="blue") +
  geom_line(data=dat, aes(x=gen, y=fZ, group=seedF), color="red") +
  geom_point(data=dat[gen==10], aes(x=gen, y=fZ), color="red") +
  geom_point(data=dat[gen==11], aes(x=gen, y=fZ), color="green") +
  geom_point(data=dat[gen==12], aes(x=gen, y=fZ), color="red")


  dat.ag <- dat[,list(delta.p=fracA[-1] - fracA[-length(neut)], p=fracA[-length(neut)], gen=gen[-length(gen)]), list(seedF)]
  ggplot(data=dat.ag[gen>200], aes(x=p, y=delta.p)) + geom_point()





  parameters <- c(c=.0000, mZZ=.5, mZW=.5)
  state <- c(U=1, V=1, W=97, X=1)/100

  gyno<-function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      dU <- c*V*mZZ/(U+V+W+X) - U
      dV <- (c*V*(1-mZZ) + (1-c)*U*V/(U+V+W+X) + (1-c)*(W*X/3)/(U+V+W+X) + (1-c)*(1/2)*(U*W + V*X)/(U+V+W+X)) - V
      dW <- c*X*mZW/(U+V+W+X) - W
      dX <- (c*X*(1-mZW)/(U+V+W+X) + (1-c)*(W*X*2/3)/(U+V+W+X) + (1-c)*(1/2)*(U*W + V*X)/(U+V+W+X)) -X

      list(c(dU, dV, dW, dX))
    })
  }

  times <- seq(0, 200, by = .1)
  library(deSolve)
  library(data.table)
  out <- as.data.table(lsoda(y=state, times = times, func = gyno, parms = parameters))

  out[,ZZ:=U+V]
  out[,ZW:=W+X]
  out[,fZ:=(ZZ+.5*ZW)/(ZZ+ZW)]

  plot(X~time, out)


  U.prime <- function(c, mZZ, mZW, U, V, W, X) V*mZZ
  V.prime <- function(c, mZZ, mZW, U, V, W, X) (V*(1-mZZ) + (1-c)^2*(U*V + (W*X)/3 + (U*X)/2 + (V*W)/2))
  W.prime <- function(c, mZZ, mZW, U, V, W, X) X*mZW
  X.prime <- function(c, mZZ, mZW, U, V, W, X) (X*(1-mZW) + (1-c)^2*((W*X*2)/3 + (U*X)/2 + (V*W)/2))

  fun <- function(c.i, mZZ.i, mZW.i) {
    mat <- matrix(NA, nrow=200, ncol=5)
    mat[,1] <- c(1:dim(mat)[1])
    mat[1,2:5] <- c(.25, .25, .25, .25)
    #c.i <- 0.01
    #mZZ.i <- .05
    #mZW.i <- .15

    for(i in 2:dim(mat)[1]) {
      mat[i,2:5] <- c(U.prime(c=c.i, mZZ=mZZ.i, mZW=mZW.i, U=mat[i-1, 2], V=mat[i-1, 3], W=mat[i-1, 4], X=mat[i-1, 5]),
                      V.prime(c=c.i, mZZ=mZZ.i, mZW=mZW.i, U=mat[i-1, 2], V=mat[i-1, 3], W=mat[i-1, 4], X=mat[i-1, 5]),
                      W.prime(c=c.i, mZZ=mZZ.i, mZW=mZW.i, U=mat[i-1, 2], V=mat[i-1, 3], W=mat[i-1, 4], X=mat[i-1, 5]),
                      X.prime(c=c.i, mZZ=mZZ.i, mZW=mZW.i, U=mat[i-1, 2], V=mat[i-1, 3], W=mat[i-1, 4], X=mat[i-1, 5]))
      mat[i,2:5] <- mat[i,2:5]/sum(mat[i,2:5])
    }

    mat <- as.data.table(mat)
    setnames(mat, names(mat), c("time", "U", "V", "W", "X"))
    mat[,ZZ:=U+V]
    mat[,ZW:=W+X]
    mat[,fZ:=(ZZ+.5*ZW)]
    mat[,male:=U+W]

    data.table(fZ=mat$fZ[dim(mat)[1]], fM=mat$U[dim(mat)[1]] + mat$W[dim(mat)[1]],
              c=c.i, mZZ=mZZ.i, mZW=mZW.i)
  }

  o <- foreach(ci=c(.01, .51, .991), .combine="rbind")%do%{
    foreach(mZZi=seq(from=0, to=1, by=.05), .combine="rbind")%do%{
      foreach(mZWi=seq(from=0, to=1, by=.05), .combine="rbind")%do%{
        fun(c.i=ci, mZZ.i=mZZi, mZW.i=mZWi)
      }
    }
  }


  x11()
  a.plot <- ggplot(data=o, aes(x=mZW, y=mZZ, fill=fM)) + geom_raster() + facet_wrap(~c) +
  scale_fill_viridis(limits=c(0,1))

  plot_grid(eq.freq.plot, a.plot)





  )
  }












  delta.fun <- function(p) {
    mean(abs(p[-length(p)] - p[-1]))
  }

  var.fun <- function(p) {
    var((p[-1] - p[length(p)]))
  }

  dat.ag <- dat[gen>400, list(delta=c(delta.fun(fZ), delta.fun(neut)),
                                gr=rep(c("fZ", "fN"))),
                        list(seedF)]

boxplot(delta~gr, dat.ag)



  var.fun <- function(p) {
    ((p[-1] - p[length(p)]))
  }

  dat.ag <- dat[gen>200, list(delta.p=c(var.fun(fracA), var.fun(neut)),
                              group=rep(c("fracA", "neut"), each=length(fracA)-1),
                              gen=c(gen[-1], gen[-1])),
                        list(seedF)]

  ggplot(data=dat.ag, aes(x=group, y=delta.p)) + geom_boxplot()



  dat.ag <- dat[n>210, list(delta.p=mean(fracA[-length(fracA)] - fracA[-1])), list(p=round(fracA, 2), seedF)]
  dat.ag <- dat[n>210, list(delta.p=mean(neut[-length(neut)] - neut[-1])), list(p=round(neut, 4), seedF)]






  delta.plot <- ggplot(data=dat.ag, aes(x=type, y=delta.p)) + geom_boxplot()
  var.plot <- ggplot(data=dat.ag, aes(x=type, y=v)) + geom_boxplot()
  plot_grid(delta.plot, var.plot)

  a.acf <- acf(dat[gen>200][seedF==1]$fracA, lag.max=400)
  neut.acf <- acf(dat[gen>200][seedF==1]$neut, lag.max=400)

verifyMarkovProperty(dat[gen>200][seedF==1]$neut)

temp <- fft(dat[gen>200][seedF==1]$neut












  dt.ag <- dt[gen>75, list(CF=median(frac_Clonal_Female), CM=median(frac_Clonal_Male),
                          SF=median(frac_Sexual_Female), SM=median(frac_Sexual_Male)),
                    list(N=param_N, CR=param_CR, AFR=param_AFR, BFR=param_BFR)]

  ggplot(data=dt.ag, aes(x=CR, y=CF+CM)) + geom_point() + geom_abline(intercept = 0, slope = 1)




  dt.ag[,maleRate_ZW:=1-AFR]
  dt.ag[,maleRate_ZZ:=1-BFR]

  eq.freq.plot <- ggplot(data=dt.ag[CR<.05 | CR>.99 | CR==.51], aes(x=maleRate_ZW, y=maleRate_ZZ, fill=(fracB))) +
  geom_raster() + scale_fill_viridis(limits=c(0,1)) + facet_wrap(~CR) + geom_abline(intercept=0, slope=1)  +
  geom_point(data= dt.ag[CR<.05 | CR>.99 | CR==.51][maleRate_ZZ==.86 & maleRate_ZW==.96], aes(x=maleRate_ZW, y=maleRate_ZZ), color="red") +
  geom_point(data= dt.ag[CR<.05 | CR>.99 | CR==.51][fracA>=(2*N-1)/(2*N) | fracA<=(1/(2*N))], aes(x=maleRate_ZW, y=maleRate_ZZ), color="blue", size=.5)

  sex.freq.plot <- ggplot(data=dt.ag[CR<.05 | CR>.99 | CR==.51], aes(x=maleRate_ZW, y=maleRate_ZZ, fill=(CM+SM))) +
  geom_raster() + scale_fill_viridis(limits=c(0,1)) + facet_wrap(~CR) + geom_abline(intercept=0, slope=1)  +
  geom_point(data= dt.ag[CR<.05 | CR>.99 | CR==.51][maleRate_ZZ==.86 & maleRate_ZW==.96], aes(x=maleRate_ZW, y=maleRate_ZZ), color="red") +
  geom_point(data= dt.ag[CR<.05 | CR>.99 | CR==.51][fracA>=(2*N-1)/(2*N) | fracA<=(1/(2*N))], aes(x=maleRate_ZW, y=maleRate_ZZ), color="blue", size=.5)




  eq.freq.plot



  dynamics.plot <- ggplot(data=dt[param_CR<.05 | param_CR>.99 | param_CR==.51][param_BFR==.86 & param_AFR==.96],
                          aes(x=gen, y=fracB, color=as.factor(fracB==0 | fracB==1))) +
  geom_point() + facet_wrap(~param_CR)





  ggplot(data=dt.ag[BFR==.86 & AFR==.96], aes(x=CR, y=fracB)) + geom_line() + geom_point() +
  geom_hline(yintercept=c(0,1))

  dt.ag[,x:=AFR>=BFR]
  ggplot(data=dt.ag, aes(x=(CR), y=fracB, group=paste(AFR, BFR, sep="_"), color=abs(AFR-BFR))) + geom_line() +
  scale_color_viridis(limits=c(0,1)) + facet_wrap(.~x) +
  geom_line(data=dt.ag[BFR==.86 & AFR==.96], aes(x=CR, y=fracB), color="red")


  plot_grid(eq.freq.plot, dynamics.plot, nrow=2)


  +
  theme(legend.position = "none")





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
