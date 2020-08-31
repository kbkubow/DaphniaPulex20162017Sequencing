### Figure 5: QTL analsys
1. `pooledMF/` performs the G' test on Pool-Seq data. Also has scripts to estimate haplotype frequencies in Pool-Seq data
2. `rQTL/` generates phased A,C, and F1s and preps data for rQTL, runs rQTL

### Figure 6: history
1. `history/` some analysis of the history of these loci

### Supplemental Figure X1:
1. `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/pooledMF_r_hybrids_questionmark.R` are MF pools really F1 pools?

### Supplemental Figure X2:
1. `DaphniaPulex20162017Sequencing/AlanAnalysis/LEA/Readme.Md` Introgression analysis

### Extraneous bit 1:
1. `DaphniaPulex20162017Sequencing/AlanAnalysis/slimuations` Alan's older slimulation & analytic model for dynamics of this locus




















1. AB_phasing/ uses whatshapp to trio data to phase A & B haplotypes
2. popPhasing/ uses shapeit to population/readbacked + priors from AB_phasing
3. pooledMF/ does sliding window test (odds_ratio_SWR.R is current one)
4. HWE_scripts/ examines HW deviations, ternary plots, clumping, etc. From ESEB
5. pulicaria/
6. slimulations/ performs forward simulations & analytical models

7. pooledMF/QTL*.R scripts perform the G' test.
8. rQTL/




need to do:

move plotting function of odds_ratio_SWR.R script into new folder/file which summarizes different bits
