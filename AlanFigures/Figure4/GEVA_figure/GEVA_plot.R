library(tidyverse)
library(magrittr)
library(reshape2)
library(forcats)
library(ggrepel)

## GEVA plotting part 2
rm(list=ls())

############################
## Load guide file
############################

load("/project/berglandlab/alan/gprime_peaks.replicates.250K.05.Rdata")
peaks_mrk = peaks
peaks_mrk %<>% 
  mutate(SNP_ID = paste(CHROM,posMaxGprime, sep = "_")) %>%
  mutate(QTL = paste("QTL", peaks_mrk$start,peaks_mrk$end, sep ="_" ) )

# Load Data
Killwood.output <- read.delim2("/scratch/yey2sn/DAPHNIA_CASE_STUDIES/GEVA_figure/Killwood_samples.txt.TMRCA.txt", sep = "\ " ,quote="\"", comment.char="", header = F)
Dorset.output <- read.delim2("/scratch/yey2sn/DAPHNIA_CASE_STUDIES/GEVA_figure/Dorset_samples.txt.TMRCA.txt", sep = "\ " ,quote="\"", comment.char="", header = F)
England.output <- read.delim2("/scratch/yey2sn/DAPHNIA_CASE_STUDIES/GEVA_figure/England_samples.txt.TMRCA.txt", sep = "\ " ,quote="\"", comment.char="", header = F)

TMRCA.output = rbind(Killwood.output,
                     Dorset.output,
                     England.output)



names(TMRCA.output) = c("CHROM","POS","ID","ClockModel","Filter","N_Concordant", "N_Discordant", "PMean","PMode","PMedian","Set")

TMRCA.output$Set = gsub("../|.txt|_samples","" , TMRCA.output$Set)


TMRCA.output %<>% mutate(SNP_ID = paste(CHROM,POS, sep = "_"))


TMRCA.output[,c("POS", "PMean","PMode","PMedian")] = lapply(TMRCA.output[,c("POS", "PMean","PMode","PMedian")] , as.numeric)

lapply(TMRCA.output , class)

############################
# Get summary stats
############################

TMRCA.output %>% 
  group_by(Set) %>%
  summarise(Mean = mean(PMean),
  SD = sd(PMean))

TMRCA.output[which(TMRCA.output$SNP_ID %in% peaks_mrk$SNP_ID),] %>% 
  group_by(Set) %>%
  summarise(Mean = mean(PMean),
            SD = sd(PMean))

TMRCA.output[which(TMRCA.output$SNP_ID %in% peaks_mrk$SNP_ID),c("SNP_ID","PMean")] 

TMRCA.output %<>%
  mutate(organize_id = ifelse(.$Set == "England", 3,
                              ifelse(Set == "Killwood", 1, 2)
                               ))



TMRCA.output %<>% 
  separate(CHROM, into = c("O_scaff","O_scaff_N", "H_scaff", "H_scaff_N"), sep = "_" , remove = F)

TMRCA.output %<>%
  mutate(QTLname = paste( paste("S",O_scaff_N, sep=""),POS, sep = "_") )

############################
# PLot distribution
############################
ggplot() + 
  geom_violin(data=TMRCA.output, aes(x=fct_reorder(Set,organize_id), y=log10(PMean))) +
  geom_boxplot(data=TMRCA.output, aes(x=fct_reorder(Set,organize_id), y=log10(PMean)),
               width=0.1, outlier.shape = NA) +
  geom_point(data=TMRCA.output[which(TMRCA.output$SNP_ID %in% peaks_mrk$SNP_ID),], 
             aes(x=fct_reorder(Set,organize_id), y=log10(PMean)), 
             color = "red",
             shape = 18,
             size = 3) +
  geom_label_repel(data=TMRCA.output[which(TMRCA.output$SNP_ID %in% peaks_mrk$SNP_ID & TMRCA.output$Set == "Killwood" ),],
    aes(x=fct_reorder(Set,organize_id), y=log10(PMean), 
            label = QTLname),
    size = 1.8,
    box.padding = 1,
    segment.color = "red") +
  ylab("Log10(TMRCA)") +
  xlab("Population Set") -> GEVA_ages

ggsave(GEVA_ages, width = 5, height = 4, file = "GEVA_ages.pdf")
  
  
############################
# PLot time manhattan
############################

ggplot() +
  geom_vline(data = peaks_mrk,
             aes(xintercept=posMaxGprime), color = "red") +
  geom_line(data =TMRCA.output,
             aes(x=POS,
                 y=PMean/100000,
                 )) +
  ylab("TMRCA (x 100,000 years)") +
  theme_classic() + 
  facet_wrap(Set~CHROM, ncol  = 2, scales = "free") -> Time_plot

ggsave(Time_plot, width = 8, height = 5, file = "Time_plot.pdf")

############################
# Plot time as distribution
############################

Peaks_TMRCA = TMRCA.output %>%
  mutate(SNP_ID = paste(CHROM,POS, sep = "_")) %>%
  .[which(.$SNP_ID %in% peaks_mrk$SNP_ID),]

Peaks_TMRCA[,c("CHROM","POS","PMean"),]

ggplot() +
  geom_histogram(data =TMRCA.output,
            aes(
                PMean/100000,
            )) +
  geom_vline(data = Peaks_TMRCA,
             aes(xintercept=PMean/100000), color = "red") +
  xlab("TMRCA (x 100,000 years)") +
  theme_classic() + 
  facet_wrap(Set~CHROM, ncol  = 2, scales = "free") -> Time_dist

ggsave(Time_dist, width = 8, height = 5, file = "Time_dist.pdf")


############################
#Plot Time Zoomed on QTLs  #
############################

tmrca_zoom = list()

for(i in 1:dim(peaks_mrk)[1]){
  
  TMRCA.output %>%
    .[which(.$CHROM == peaks_mrk$CHROM[i]),] -> tmp
  
  tmp %>%
    .[which(.$POS > peaks_mrk$start[i]  & .$POS  < peaks_mrk$end[i] ),] %>% 
    mutate(QTL = paste("QTL", peaks_mrk$start[i],peaks_mrk$end[i], sep ="_" ) ) -> QTL_slice
  
  tmrca_zoom[[i]] = QTL_slice
}

tmrca_zoom_df = do.call(rbind, tmrca_zoom)

## PLOT ##

ggplot() +
  geom_vline(data = peaks_mrk,
             aes(xintercept=posMaxGprime), color = "red") +
  geom_line(data =tmrca_zoom_df,
            aes(x=POS,
                y=PMean/100000,
            )) +
  ylab("TMRCA (x 100,000 years)") +
  theme_classic() + 
  facet_wrap(~QTL, ncol  = 2, scales = "free") -> Time_plot_zoom

ggsave(Time_plot_zoom, width = 6, height = 6, file = "Time_plot_zoom.pdf")


############################
#Sliding Window approach   #
############################

CHRS = unique(TMRCA.output$CHROM)
rolling_TMRCAS_list = list()

for(i in 1:length(CHRS) ) {
  
  tmp = TMRCA.output[which(TMRCA.output$CHROM == CHRS[i] ),]
  
  First_pos = min(tmp$POS)
  Last_pos = max(tmp$POS)
  Complete_window = seq(from = First_pos, to=Last_pos, by =1 )
 
   Artificial_chromosome = data.frame(POS=Complete_window) %>% 
    full_join(., tmp[,c("POS", "PMean")])
   
   Rolling_Posterior = slide( 
                 .x= Artificial_chromosome$PMean, 
                 .f = function(x) mean(x, na.rm = T), 
                 .after = 250000,
                 .step = 50000, 
                 .complete = FALSE)
 
   Rolling_Position = slide( 
                 .x= Artificial_chromosome$POS, 
                 .f = function(x) mean(x, na.rm = T), 
                 .after = 250000,
                 .step = 50000, 
                 .complete = FALSE)
   
   
   data.frame(CHROM=CHRS[i] , 
              POS = do.call(rbind,Rolling_Position), 
              ROLLING_TMRCA= do.call(rbind,Rolling_Posterior)
              ) -> tmp_out
  
   rolling_TMRCAS_list[[i]] = tmp_out

} # i loop

rolling_TMRCAS_df = do.call(rbind,rolling_TMRCAS_list)


############################
# Plot Sliding Window      #
############################

ggplot() +
  geom_vline(data = peaks_mrk,
             aes(xintercept=posMaxGprime/100000), color = "red") +
  geom_line(data =rolling_TMRCAS_df,
            aes(x=POS/100000,
                y=ROLLING_TMRCA/100000
            )) +
  ylab("Mean TMRCA (x 100,000 years)") +
  xlab("Position (x 100,000 bp)") +
  theme_bw() + 
  facet_wrap(~CHROM, ncol  = 3, scales = "free") -> Sliding_Time_plot

ggsave(Sliding_Time_plot, width = 8, height = 5, file = "Sliding_Time_plot.pdf")
