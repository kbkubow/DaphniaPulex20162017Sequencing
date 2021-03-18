### libraries
  library(data.table)
  library(ggplot2)
  library(viridis)
  library(tidyverse)
  library(ggrepel)
  library(cowplot)
  library(patchwork)

### set working directory
  setwd("/Users/alanbergland/Documents/GitHub")

############
### data ###
############

### load IBS data
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure2/ibs.long.Rdata") ### loads `ibs.long` object

### load diversity estimates (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure2/Theta_pi_clone.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure2/genome_theta.Rdata")

### load ROH information (made by )
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure2/rohan_summary.Rdata")


#######################
### plot components ###
#######################

### IBS
  h.just <- .25
  v.just <- .25
  l.size <- 1.5

  corrmatrix <- ggplot(data=ibs.long, aes(scid.a, scid.b, fill=IBS)) +
                geom_raster() +
                scale_fill_viridis(option="D") +
                theme(axis.title=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank())

### Theta-pi
  ## among ponds
    theta_ponds <-
      ggplot(out[clone %in% c(rand.sc, oo.sc)][data %in% c("theta_inc_ROH")],
           aes(x=year, y=global)) +
      geom_boxplot() +
      facet_wrap(~factor(pond, levels=c("DCat", "D8", "DBunk"))) +
      xlab("Year") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
            axis.text.y = element_text(size=12),
            axis.title = element_text(size=14)) +
      scale_y_continuous(breaks=c(.002, .0025, .003, .0035), labels=c(2, 2.5, 3, 3.5), limits=c(0.0019, 0.00375)) +
      ylab(expression(paste(theta[pi], " (x", 10^-3, ")", sep="")))


  # AxC crosses
    theta_AC <-
    ggplot() +
    geom_point(data=out[!clone %in% (cross$clone)][SC %in% c("A", "C")][data %in% c("theta_inc_ROH")][pond=="D8"],
               aes(x=SC, y=global), size=2) +
    geom_point(data=out[clone %in% cross[cross=="F1"]$clone][!SC=="AL"][data %in% c("theta_inc_ROH")][pond=="D8"],
               aes(x="F1-AxC", y=global), size=2) +
    geom_point(data=out[clone %in% cross[cross=="F2"]$clone][data %in% c("theta_inc_ROH")][pond=="D8"],
               aes(x="F2-AxC", y=global), size=2) +
    theme_bw() +
    facet_wrap(~pond) +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=14),
          legend.position = "none") +
    scale_y_continuous(breaks=c(.002, .0025, .003, .0035), labels=c(2, 2.5, 3, 3.5), limits=c(0.0019, 0.00375)) +
    ylab(expression(paste(theta[pi], " (x", 10^-3, ")", sep=""))) +
    xlab("")

# Selfed vs offspring (Supp Fig)
  theta_selfed <-
    ggplot() +
    geom_point(data=out[clone %in% selfing.par][data %in% c("theta_inc_ROH")],
               aes(x=as.factor(c("1", "1", "1", "1")),
                   y=global), size=2) +
    geom_point(data=out[SC %in% c("B", "H", "C", "W")][data %in% c("theta_inc_ROH")],
               aes(x=as.factor(SC), y=global), size=2) +
    labs(x="", y="") +
    scale_x_discrete(labels=c("B"="Offspring", "C"="Offspring",
                              "H"="Offspring", "W"="Offspring",
                              "1"="Parent")) +
    facet_grid(~SC.share, scales = "free_x") +
    guides(color=FALSE) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size=12, angle = 45, hjust=1, ),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=14)) +
    scale_y_continuous(breaks=c(.002, .0025, .003, .0035), labels=c(2, 2.5, 3, 3.5), limits=c(0.0019, 0.00375)) +
    ylab(expression(paste(theta[pi], " (x", 10^-3, ")", sep="")))



### ROHs
  rl.ag <- rl[ROH_LENGTH>=100000, list(nroh=.N, sroh=sum(ROH_LENGTH), mroh=mean(ROH_LENGTH)), list(clone, SC.uniq, population)]
  rl.ag.ag <- rl.ag[, list(nroh=mean(nroh), sroh=mean(sroh), mroh=mean(mroh)), list(SC.uniq, population)]


  rl.ag.ag[SC.uniq%in%c("A", "C") & population=="D8",lab:=SC.uniq]

  rl.ag.ag[is.na(lab), lab:=""]
  rl.ag.ag[population=="Dcat", population:="DCat"]

  m.ag <- m[,list(perc_roh=mean(perc_roh), len_roh=mean(len_roh)), list(SC.uniq, population)]
  m.ag[SC.uniq%in%c("A", "C") & population=="D8",lab:=SC.uniq]

  m.ag[is.na(lab), lab:=""]

  t1 <- lm(nroh~sroh, rl.ag.ag[sroh<2.5e6])


  summary(lm(sroh~population, rl.ag.ag[population%in%c("D8", "DBunk", "DCat")]))

  roh_plot <-
  ggplot(data=rl.ag.ag[population%in%c("D8", "DBunk", "DCat")],
        aes(x=sroh, y=nroh, label=lab)) +
  geom_abline(aes(slope=coef(t1)[2], intercept=coef(t1)[1]))  +
  geom_point() +
  geom_label_repel(box.padding=1.5) +
  facet_wrap(~factor(population, levels=c("DCat", "D8", "DBunk"))) +
  theme_bw() +
  ylab("Number ROH") +
  xlab(expression(paste("Summed length ROH (x", 10^6, ")", sep=""))) +
  theme(legend.position = "none",
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14)) +
  scale_x_continuous(breaks=c(0, 2.5e6, 5e6, 7.5e6), labels=c(0, 2.5, 5, 7.5))




### pedigree
  pedigree_plot <-
  ggplot() +
  draw_image("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure2/Pedigree.png") +
  theme(axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())



#################
### MEGA-plot ###
#################

layout <- "
AABBC
AAEED
FFFFF
FFFFF
FFFFF"


layout <- "
AABBEE
AACDEE
FFFFFF
FFFFFF
FFFFFF"


layout <- "
AAFFFF
AAFFFF
BBBBEE
CCDDEE
"




mega <-
corrmatrix + theta_ponds + theta_AC + theta_selfed + roh_plot + pedigree_plot +
theme(strip.text=element_text(size=12), plot.tag = element_text(size=16)) +
plot_layout(design = layout) +
plot_annotation(tag_levels = 'A')


ggsave(mega, file="~/mega_diversity_v2.pdf")
