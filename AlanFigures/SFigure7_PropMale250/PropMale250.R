#!/usr/bin/env Rscript

### libraries
  library(data.table)
  library(ggplot2)

### Load in data files
  grace250spring <- fread("/Users/kbkubow/Box Sync/UVA/Daphnia/A_B_1liters/Grace250Sp2020.csv")
  grace250spring$SC <- ifelse(grace250spring$Clone=="D8179" | grace250spring$Clone=="D8349", "A", "C")
  grace250spring$momsplusslide <- grace250spring$Males_4 + grace250spring$Moms_4 + grace250spring$FemalesSlide_4
  grace250spring$propmale <- grace250spring$Males_4/grace250spring$momsplusslide

### Graph
  grace250spring[,Clone := factor(Clone, levels=c("D8179", "D8349", "D8222", "D8515"))]
  Grace250 <- ggplot(data=grace250spring, aes(x=Clone, y=propmale, color=SC, group=Clone)) +
    geom_boxplot() + ylim(0,0.45) + ylab("Prop Male") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
