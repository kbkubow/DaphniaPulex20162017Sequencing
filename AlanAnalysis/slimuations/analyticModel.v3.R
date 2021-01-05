### this model incorporates seasonal death of clonal individuals

#### libraries
  library(data.table)
  library(foreach)
  library(cowplot)
  library(viridis)
  library(ggplot2)
  #library(deSolve)
  library(doMC)
  registerDoMC(4)

### function
  fun <- function(c, mZZ, mZW, mWW, x, cs) {
    ### c is the probabilty of clonal reproduction
    ### mZZ, mZW, mWW are the male production rates of these genotypes
    ### x is the vector of genotype frequencies at time t
    ### return: the vector of genotype frequencies at time t+1

    U <- x[1]
    V <- x[2]
    W <- x[3]
    X <- x[4]
    Y <- x[5]
    Z <- x[6]

    U.prime <- c*mZZ*V*cs   ### Male ZZ
    V.prime <- (c*(1-mZZ)*V*cs + 2*(1-c)*(U*V + W*V/2 + U*X/2 + W*X/4))   ### Female ZZ

    W.prime <- c*(mZW)*X*cs    ### Male ZW
    X.prime <- c*(1-mZW)*X*cs + 2*(1-c)*(W*X/2 + U*X/2 + W*V/2 + W*Z/2 + Y*X/2)    ### Female ZW

    Y.prime <- c*(mWW)*Z*cs    ### Male WW
    Z.prime <- (c*(1-mWW)*Z*cs + 2*(1-c)*(Y*Z + W*X/4 + W*Z/2 + Y*X/2))           ### Female WW

    c(U.prime, V.prime, W.prime, X.prime, Y.prime, Z.prime)/(sum(c(U.prime, V.prime, W.prime, X.prime, Y.prime, Z.prime)))
  }

  sim <- function(c, mZZ, mWW, d, nGens, init, clonal_survival) {
    #c=.5
    #mZZ=.5
    #mWW=.5
    #nGens=200
    #d=.5
    #lethal=0
    #init <- rep(1,6)/6
    #clonal_survival=c(rep(1, 9), .1)

    mZW <- mZZ + d*(mWW - mZZ)

    mat <- matrix(NA, nrow=nGens, ncol=6)
    mat[1,] <- init

    for(i in 2:nGens) {
      mat[i,] <- fun(c=c, mZZ=mZZ, mZW=mZW, mWW=mWW, x=mat[i-1,], cs=clonal_survival[i%%10+1])
    }

    mat <- cbind(mat,
                (mat[,1]+mat[,2]) + .5*(mat[,3] + mat[,4]))

    mat.dt <- data.table(freq=expand.grid(mat)[,1], gen=rep(c(1:nGens), 7),
                        class=rep(c("m.ZZ", "f.ZZ", "m.ZW", "f.ZW", "m.WW", "f.WW", "fZ"), each=nGens),
                        c=c, d=d, mZZ=mZZ, mZW=mZW, mWW=mWW,
                        cs=clonal_survival[rep(c(1:nGens), 7)%%10+1])

    #ggplot(mat.dt, aes(x=gen, y=freq, group=class, color=class)) + geom_line() + facet_grid(~I(class=="fZ")) + ylim(0,1)

    return(mat.dt)
  }

### simulation
  o <- foreach(c.i=c(0.01, .5, .991), .combine="rbind")%dopar%{
    print(c.i)
    foreach(mZZ.i=seq(from=.01, to=.99, by=.05), .combine="rbind")%do%{
      foreach(mWW.i=seq(from=.01, to=.99, by=.05), .combine="rbind")%do%{
        foreach(d.i=c(.5), .combine="rbind")%do%{
          #c.i=.5; mZZ.i=.5; mWW.i=.5; d.i=.5; l.i=0
            message(paste(c.i, mZZ.i, mWW.i, d.i, sep=" / "))
            o.tmp <- sim(c=c.i, mZZ=mZZ.i, mWW=mWW.i, d=d.i,
                    clonal_survival=c(rep(1, 9), 0),
                    nGens=1000,
                    init=rep(1,6)/6)
            #o.tmp[gen>500, list(freq=mean(freq), var.freq=var(freq)), list(c, d, mZZ, mZW, mWW, class)]
            o.tmp
        }
      }
    }
  }
  o.ag <-o[gen>500, list(freq=mean(freq), var.freq=var(freq)), list(c, d, mZZ, mZW, mWW, class)]


  ggplot(data=o[class=="fZ"], aes(x=gen, y=freq, group=interaction(mZZ, mWW))) + geom_line() + facet_wrap(~c)

  o[gen>900][c==.5][class=="fZ"][freq>.6 & freq<.66]
  o.ag[class=="fZ"][d==.5 & mZZ==.01 & mZW==.26 & mWW==.51]

  ggplot(data=o[class=="fZ"][d==.5 & mZZ==.01 & mZW==.26 & mWW==.51],
          aes(x=gen, y=freq, group=interaction(mZZ, mWW))) + geom_line() + facet_wrap(~c)
  o.ag[class=="fZ"][d==.5 & mZZ==.01 & mZW==.26 & mWW==.51]

  ggplot(data=o.ag[class=="fZ"], aes(x=mWW, y=mZZ, fill=freq)) +
  geom_raster() + facet_grid(d~c) + scale_fill_viridis(limits=c(-.05,1.05)) +
  geom_abline(intercept = 0, slope = 1) + ggtitle("Overwintering death of clonally produced individuals every 10 gens")

  ggplot(data=o.ag[class=="fZ"], aes(x=mWW, y=mZZ, fill=freq/sqrt(var.freq))) +
  geom_raster() + facet_grid(d~c) + scale_fill_viridis() +
  geom_abline(intercept = 0, slope = 1) + ggtitle("Overwintering Death every 10 gens")



fixation <- -0.5

ggplot(data=o[class=="fZ"][freq>fixation & freq<(1-fixation)], aes(x=mWW, y=mZZ, fill=freq)) +
geom_raster() + facet_grid(d~c) + scale_fill_viridis(limits=c(-.05,1.05)) +
geom_abline(intercept = 0, slope = 1) + ggtitle("Overwintering Death every 10 gens")









  loadSlim <- function(fn) {
    dat <- fread(fn)
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
    dat
  }

  loadSlim.obj <- function(fn) {
    load(fn)
    dat <- dt
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
    dat
  }

  slim_analytic <- function(dat) {
    dat.temp <- dat[gen==500]
    setkey(dat.temp, N, CR, AMR, BMR)
    dat.temp <- dat.temp[!duplicated(dat.temp)]

    sim.dat <- foreach(i=1:dim(dat.temp)[1], .combine="rbind")%do%{

      sim.dat  <- sim(c=dat.temp[i]$CR,
                      mZW=dat.temp[i]$AMR,
                      mZZ=dat.temp[i]$BMR,
                      nGens=490,
                      init=c(0.01, .01, 0.01, 0.97))[class=="fZ",c("freq", "gen"), with=F]
      sim.dat[,gen:=gen+10]
      setnames(sim.dat, names(sim.dat), c("fZ", "gen"))
      sim.dat[, seedF:= -1*(dat.temp[i]$seedF)]
      sim.dat[, CR:=dat.temp[i]$CR]
      sim.dat[, AMR:=dat.temp[i]$AMR]
      sim.dat[, BMR:=dat.temp[i]$BMR]
      sim.dat
    }

    dat <- rbind(sim.dat, dat, fill=T)

    dat
  }




  fixation <- 0.01

  lethal_1 <- ggplot(data=o[class=="fZ"][lethal=="0;1"][freq>fixation & freq<(1-fixation)], aes(x=mWW, y=mZZ, fill=freq)) +
  geom_raster() + facet_grid(d~c ) + scale_fill_viridis(limits=c(-.05,1.05)) +
  geom_abline(intercept = 0, slope = 1) + ggtitle("WW lethal")

  lethal_01 <- ggplot(data=o[class=="fZ"][lethal=="0.01;0.01"][freq>fixation & freq<(1-fixation)], aes(x=mWW, y=mZZ, fill=freq)) + geom_raster() + facet_grid(d~c ) + scale_fill_viridis(limits=c(-.05,1.05)) +
  geom_abline(intercept = 0, slope = 1) + ggtitle("1% inbreeding depression")


  lethal_05 <- ggplot(data=o[class=="fZ"][lethal=="0.05;0.05"][freq>fixation & freq<(1-fixation)], aes(x=mWW, y=mZZ, fill=freq)) + geom_raster() + facet_grid(d~c ) + scale_fill_viridis(limits=c(-.05,1.05)) +
  geom_abline(intercept = 0, slope = 1) + ggtitle("5% inbreeding depression")

  lethal_0 <- ggplot(data=o[class=="fZ"][lethal=="0;0"][freq>fixation & freq<(1-fixation)], aes(x=mWW, y=mZZ, fill=freq)) + geom_raster() + facet_grid(d~c ) + scale_fill_viridis(limits=c(-.05,1.05)) +
  geom_abline(intercept = 0, slope = 1) + ggtitle("No inbreeding depression")

  lethal_0 + lethal_01 + lethal_05 + lethal_1  + plot_layout(guides = 'collect')



### contrast simulation to analytical model
  #o <- slim_analytic(loadSlim("~/slim_out.delim"))
  o <- slim_analytic(loadSlim.obj("~/dt.Rdata"))

  o.ag.exp <-  o[seedF<0, list(fZ.exp=mean(fZ)),
                                    list(gen, CR, AMR, BMR)]

  o.ag.obs <-  o[seedF>0, list(fZ.obs=mean(fZ)),
                                    list(gen, CR, AMR, BMR)]

  setkey(o.ag.exp, gen, CR, AMR, BMR)
  setkey(o.ag.obs, gen, CR, AMR, BMR)
  o.ag <- merge(o.ag.exp, o.ag.obs)

  o.ag.ag <- o.ag[gen>300, list(fZ.exp=mean(fZ.exp), fZ.obs=mean(fZ.obs)), list(CR, AMR, BMR)]

  #o.ag.ag <- o.ag[gen>450, list(fZ.exp=mean(fZ.exp), fZ.delta=fZ.obs[-1] - fZ.obs[-length(fZ.obs)], fZ=fZ.obs[-length(fZ.obs)]), list(CR, AMR, BMR)]


  ind.plot <- ggplot() +
              geom_line(data=o[gen>200], aes(x=gen, y=neut, group=seedF), color="blue") +
              geom_line(data=o[seedF>0], aes(x=gen, y=fZ, group=seedF), color="red") +
              geom_line(data=o[seedF<0], aes(x=gen, y=fZ, group=seedF), color="black") +
              ylab("Frequency") + xlab("Generation") + ylim(0,1)


  mu.plot <- ggplot(data=o.ag) +
            geom_line(aes(x=gen, y=fZ.exp, group=paste(CR, AMR, BMR, sep="_")), color="black") +
            geom_line(aes(x=gen, y=fZ.obs, group=paste(CR, AMR, BMR, sep="_")), color="red") +
            xlab("Generation") + ylab("") + ylim(0,1)

  pw.plot <- ggplot(data=o.ag.ag, aes(x=fZ.exp, y=I(fZ.obs))) + geom_point() + geom_abline(slope=1, intercept=0)
  #plot_grid(ind.plot, mu.plot)
  pw.plot + geom_vline(xintercept=2/3)


  o.ag.ag[,delta:=fZ.obs - fZ.exp]
  summary(lm((delta)~CR*AMR*BMR , o.ag.ag))





  load("~/dt.Rdata")
  dt[,fZW := (nZW / n)]
  dt[,fZZ := (nZZ / n)]

  dt.ag <- dt[gen>200, list(fZW=mean(fZW), fZZ=mean(fZZ), fZ=mean(fZ),
                          CF=mean(frac_Clonal_Female), CM=mean(frac_Clonal_Male),
                          SF=mean(frac_Sexual_Female), SM=mean(frac_Sexual_Male),
                          neutFixation=any(neut==0 | neut==1),
                          neutDelta=neut[gen==500] - neut[gen==201],
                          fZdelta=fZ[gen==500] - fZ[gen==201]),
                    list(N=N, CR=CR, maleRate_ZW=AMR, maleRate_ZZ=BMR)]


  o <- foreach(i=1:dim(dt.ag)[1], .combine="rbind")%do%{
      print(i)
      sim.o <- sim(c=dt.ag[i]$CR, mZZ=dt.ag[i]$maleRate_ZZ, mZW=dt.ag[i]$maleRate_ZW,
          nGens=50, init=c(.25, .25, .25, .05))[gen==50]

      tmp <- dt.ag[i]
      tmp[,exp.fZ:=sim.o[class=="fZ"]$freq]
      tmp
  }


  ggplot(data=o, aes(x=fZ, y=exp.fZ)) + geom_point() + facet_wrap(~CR) + ylim(0,1) + xlim(0,1)

  o <- sim(c=.01, mZZ=.01, mZW=.15,
          nGens=50, init=c(.25, .25, .25, .25))

  ggplot(data=o, aes(x=gen, y=freq, group=class, color=class)) + geom_line()
