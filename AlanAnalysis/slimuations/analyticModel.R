#### libraries
  library(data.table)
  library(foreach)
  library(cowplot)
  library(viridis)
  library(ggplot2)
  #library(deSolve)

### function
  fun <- function(c, mZZ, mZW, x) {
    W <- x[1]
    X <- x[2]
    Y <- x[3]
    Z <- x[4]

    W.prime <- c*mZZ*X
    X.prime <- c*(1-mZZ)*X + 2*(1-c)*(W*X + X*Y/2 + Z*Y/3 + Z*W/2)
    Y.prime <- c*(mZW)*Z
    Z.prime <- c*(1-mZW)*Z + 2*(1-c)*(X*Y/2 + 2*Z*Y/3 + Z*W/2)

    c(W.prime, X.prime, Y.prime, Z.prime)/(sum(c(W.prime, X.prime, Y.prime, Z.prime)))
  }

  sim <- function(c, mZZ, mZW, nGens, init) {
    #c=.01
    #mZZ=.15
    #mZW=.01
    #nGens=200
    #  init <- c(.25, .25, .25, .25)

    mat <- matrix(NA, nrow=nGens, ncol=4)
    mat[1,] <- init

    for(i in 2:nGens) {
      mat[i,] <- fun(c=c, mZZ=mZZ, mZW=mZW, x=mat[i-1,])
    }

    mat <- cbind(mat,
                (mat[,1]+mat[,2]) + .5*(mat[,3] + mat[,4]))

    mat.dt <- data.table(freq=expand.grid(mat)[,1], gen=rep(c(1:nGens), 5),
                        class=rep(c("m.ZZ", "f.ZZ", "m.ZW", "f.ZW", "fZ"), each=nGens), c=c, mZZ=mZZ, mZW=mZW)



    return(mat.dt)
  }

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

### simulation
  o <- foreach(c.i=c(0.01, .51, .991), .combine="rbind")%do%{
    print(c.i)
    foreach(mZZ.i=seq(from=.01, to=.99, by=.05), .combine="rbind")%do%{
      foreach(mZW.i=seq(from=.01, to=.99, by=.05), .combine="rbind")%do%{
        sim(c=c.i, mZZ=mZZ.i, mZW=mZW.i,
                nGens=50, init=c(.25, .25, .25, .05))[gen==50]
      }
    }
  }

  ggplot(data=o[class=="fZ"], aes(x=mZW, y=mZZ, fill=freq)) + geom_raster() + facet_wrap(~c) + scale_fill_viridis() +
  geom_abline(intercept = 0, slope = 1)


### contrast simulation to analytical model
  o <- slim_analytic(loadSlim("~/slim_out.delim"))

  o.ag.exp <-  o[seedF<0, list(fZ.exp=mean(fZ)),
                                    list(gen, CR, AMR, BMR)]

  o.ag.obs <-  o[seedF>0, list(fZ.obs=mean(fZ)),
                                    list(gen, CR, AMR, BMR)]

  setkey(o.ag.exp, gen, CR, AMR, BMR)
  setkey(o.ag.obs, gen, CR, AMR, BMR)
  o.ag <- merge(o.ag.exp, o.ag.obs)


  ind.plot <- ggplot() +
  geom_line(data=o[gen>200], aes(x=gen, y=neut, group=seedF), color="blue") +
  geom_line(data=o[seedF>0], aes(x=gen, y=fZ, group=seedF), color="red") +
  geom_line(data=o[seedF<0], aes(x=gen, y=fZ, group=seedF), color="black") +
  ylab("Frequency") + xlab("Generation")


  ggplot(data=o.ag) +
  geom_line(aes(x=gen, y=fZ.exp, group=paste(CR, AMR, BMR, sep="_")), color="red") +
  geom_line(aes(x=gen, y=fZ.obs, group=paste(CR, AMR, BMR, sep="_")), color="blue") +
  ylim(0,1)











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
