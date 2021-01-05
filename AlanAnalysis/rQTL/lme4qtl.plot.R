
### scp aob2x@rivanna.hpc.virginia.edu:~/lme4qtl_output.AxC.long.Rdata ~/.
  library(data.table)
  library(ggplot2)

  load("~/lme4qtl_output.AxC.long.Rdata")

  o.ag <- o[, list(p.z=mean(p.z, na.rm=T), chisq=max(chisq), p.aov=min(p.aov, na.rm=T)), list(term, perm, id, chr, pos)]
  o.ag.ag <- o.ag[,list(pr=mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id, chr, pos)]
  o.ag.ag[pr==1, pr:=1/201]

### ANOVA output
  o.ag.perm <- o.ag[perm!=0, list(p.aov=c(quantile(p.aov, .01), quantile(p.aov, .05)),
                                  q=c(.01, .05)), list(term, chr)]

  setnames(peaks, "CHROM", "chr")

  f1.plot <- ggplot() +
  geom_line(data=o.ag[perm==0], aes(x=pos, y=-log10(p.aov), color=chr)) +
  geom_vline(data=peaks, aes(xintercept=posMaxGprime, linetype=as.factor(rep))) +
  geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.aov), linetype=as.factor(q))) +
  facet_grid(term~chr) +
  theme(legend.position = "none")



### beta-p output
  o.ag.perm <- o.ag[perm!=0, list(p.z=quantile(p.z, .05)), list(term, chr)]

  ggplot() +
  geom_line(data=o.ag[perm==0], aes(x=pos, y=-log10(p.z), color=chr)) +
  geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.z))) +
  facet_grid(term~chr)+
  theme(legend.position = "none")


### normalized
  ggplot() +
  geom_line(data=o.ag.ag, aes(x=pos, y=-log10(1-pr), color=chr)) +
  facet_grid(term~chr)+
  theme(legend.position = "none")
