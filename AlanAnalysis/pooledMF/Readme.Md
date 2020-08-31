# Readme for scripts

## Phase
## Estimate haplotype frequencies using Harp
`DaphniaPulex20162017Sequencing/AlanAnalysis/AB_phasing/Readme.Md`

## G' test 
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/aseReadCounter_mf.sh` runs ASE read-counter to call allele frequencies
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLSeqR.R` runs Gprime analysis
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLPlot.R` makes composite QLT plot, relies on HARP data

## HWE at QTL regions in F1s
###



find /home/aob2x/ncbi/public/sra -type f -mtime +7 -name '*.cache' -execdir rm -- '{}' \;
