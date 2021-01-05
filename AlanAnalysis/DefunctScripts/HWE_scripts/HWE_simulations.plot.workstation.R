### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(ggtern)
  library(viridis)

### load data
  load(file="/mnt/sammas_storage/bergland-lab/Daphnia_HWE/hwe.ag.m.ag")

### DeFinetti enrichment diagraim
  hwe.ag.m.ag <- hwe.ag.m[,list(diff.mu=mean(diff), diff.sd=sd(diff), diff.pr=mean(diff>0),
                                n.obs=mean(n.obs), n.sim=mean(n.sim), sd.sim=sd(n.sim)),
                         list(fAA, fAa, faa, py)]
  save(hwe.ag.m.ag, file="~/hwe.Rdata")
  hwe.ag.m.ag[,p:=apply(cbind(pnorm(n.obs, n.sim, sd.sim, log=T, lower.tail=F), pnorm(n.obs, n.sim, sd.sim, log=t, lower.tail=T)), 1, min)]



  hwe.ag.m.ag[py=="D8.2016.2017.2018.2019"][p!=-Inf][which.min(p)]
  hwe.ag.m.ag[,t:=diff.mu/diff.sd]

  ggsave(ggplot() +
        geom_point(data=hwe.ag.m.ag[py=="D8.2016.2017.2018.2019"][order(diff.mu, decreasing=F)][n.obs>2],
                    aes(x=fAA, y=fAa, z=faa, color=sign(diff.mu)*(abs(diff.mu))), size=1.5) +
        coord_tern(expand=T) + limit_tern(T = 1.05, L = 1.05, R = 1.05) +
        scale_color_viridis() + scale_fill_viridis(),
        file="~/D8_tern.pdf")


  

  ggsave(ggplot() +
        geom_point(data=hwe.ag.m.ag[py=="DBunk.2016.2017.2018.2019"][order(diff.mu, decreasing=F)][diff.mu!=0][n.obs>2],
                    aes(x=fAA, y=fAa, z=faa, color=sign(diff.mu)*(abs(diff.mu))), size=.85) +
        coord_tern(expand=T) + limit_tern(T = 1.05, L = 1.05, R = 1.05) +
        scale_color_viridis() + scale_fill_viridis(),
        file="~/DBunk_tern.pdf")
