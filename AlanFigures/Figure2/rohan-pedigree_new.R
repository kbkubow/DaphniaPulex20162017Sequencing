# Connor Murray
# 7.15.2020

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(cowplot)

#setwd("C:/Users/Conno/Desktop/Spring2020/rohan/data/")
setwd("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure2")

# Filter metadata sample file
samps <- fread("cloneinfo.txt")
samps <- data.table(samps %>% filter(Nonindependent==0 & is.na(LowReadDepth) & Species=="pulex"))

# Summary output from ROHan
out <- readRDS("summary.all.rds")
out[pond=="Dcat"]$pond <-"DCat"

# Subsample 1 clone per SC
rand.sc <- data.table(samps %>% filter(population %in% c("D8","DBunk","DCat")) %>%
                        group_by(SC) %>% sample_n(1))$clone

# OO clones of interest
oo.sc <- data.table(samps %>% filter(population %in% c("D8","DBunk","DCat") & SC %in% "OO"))$clone

out$pond <- factor(out$pond, levels = c("DBunk", "D8", "DCat"))

# Theta pi including ROH through time and ponds (FIG 2)
a <- ggplot(out[clone %in% c(rand.sc, oo.sc)][data %in% c("theta_inc_ROH")],
       aes(x=year, y=global)) +
  geom_boxplot() +
  facet_wrap(~pond) +
  labs(x="Year", y="Theta pi") +
  ylim(0.0017, 0.004) +
  theme_classic() +
  theme(axis.text.x =element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        title = element_text(face="bold", size=15),
        strip.text = element_text(face="bold", size=12))

# Parent-offspring relationships
par.off <- data.table(rbind(read.csv("D8ParentOffspringRelationships.long.csv"),
                            read.csv("DBunkParentOffspringRelationships.long.csv")))

# A x C cross clones
cross <- data.table(read.csv("A.C.cross.csv", header = TRUE))

# AxC crosses (Supp Fig)
b <- ggplot() +
  geom_point(data=out[!clone %in% (cross$clone)][SC %in% c("A", "C")][data %in% c("theta_inc_ROH")],
             aes(x=SC, y=global, color=SC), size=3) +
  geom_point(data=out[clone %in% cross[cross=="F1"]$clone][!SC=="AL"][data %in% c("theta_inc_ROH")],
             aes(x="F1-AxC", y=global), color="purple", size=4) +
  geom_point(data=out[clone %in% cross[cross=="F2"]$clone][data %in% c("theta_inc_ROH")],
             aes(x="F2-AxC", y=global), color="orange", size=4) +
  labs(x="", y="Theta pi") +
  ylim(0.0017, 0.004) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x =element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15))

# Selfing parents
selfing.par <- c(samps[SC %in% c("poB", "poH", "poC", "poW")]$clone)

# Add column for graph
out[clone %in% selfing.par, SC.share:=str_replace(paste(SC), "po", "")]
out[!clone %in% selfing.par, SC.share:=paste(SC)]

# Selfed vs offspring (Supp Fig)
c <- ggplot() +
  geom_point(data=out[clone %in% selfing.par][data %in% c("theta_inc_ROH")],
             aes(x=as.factor(c("1", "1", "1", "1")),
                 y=global, color=c("B", "C", "H", "W")), size=4) +
  geom_point(data=out[SC %in% c("B", "H", "C", "W")][data %in% c("theta_inc_ROH")],
             aes(x=as.factor(SC), y=global, color=SC), size=4) +
  labs(x="", y="") +
  scale_x_discrete(labels=c("B"="Offspring", "C"="Offspring",
                            "H"="Offspring", "W"="Offspring",
                            "1"="Parent")) +
  facet_grid(~SC.share, scales = "free_x") +
  ylim(0.0017, 0.004) +
  guides(color=FALSE) +
  theme_bw() +
  theme(legend.position = c(0.9, 0.9),
        legend.background = element_blank(),
        axis.text.x =element_text(face="bold", size=11),
        axis.text.y = element_text(face="bold", size=11),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        strip.text = element_text(face="bold", size=12))

# Cowplot of all graphs
plot_grid(a, plot_grid(b, c, labels = c("B.", "C."),
                       rel_widths = c(1,2), align="h", axis = "tb"),
          nrow = 2, labels = c("A.", ""))



###roh plots

a <- ggplot(out[clone %in% c(rand.sc, oo.sc)][data %in% c("theta_inc_ROH", "seg_in_ROH%")],
       aes(x=year, y=global)) +
  geom_boxplot() +
  facet_grid(data~pond, scales="free_y") +
  labs(x="Year", y="Theta pi") +
  theme_classic() +
  theme(axis.text.x =element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        title = element_text(face="bold", size=15),
        strip.text = element_text(face="bold", size=12))

        b <- ggplot() +
          geom_point(data=out[!clone %in% (cross$clone)][SC %in% c("A", "C")][data %in% c("theta_inc_ROH", "seg_in_ROH%")],
                     aes(x=SC, y=global, color=SC), size=3) +
          geom_point(data=out[clone %in% cross[cross=="F1"]$clone][!SC=="AL"][data %in% c("theta_inc_ROH", "seg_in_ROH%")],
                     aes(x="F1-AxC", y=global), color="purple", size=4) +
          geom_point(data=out[clone %in% cross[cross=="F2"]$clone][data %in% c("theta_inc_ROH", "seg_in_ROH%")],
                     aes(x="F2-AxC", y=global), color="orange", size=4) +
          labs(x="", y="Theta pi") +
          theme_bw() +
          facet_grid(data~., scales="free_y") +

          theme(legend.position = "none",
                axis.text.x =element_text(face="bold", size=12),
                axis.text.y = element_text(face="bold", size=12),
                axis.title.x = element_text(face="bold", size=15),
                axis.title.y = element_text(face="bold", size=15))
