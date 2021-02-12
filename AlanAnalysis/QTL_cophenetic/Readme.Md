# Make trees using iqtree w/ bootstraps (bionj -fast option):
  ## `DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_cophenetic/analysis.makeTrees.sh`

# collect trees and write summary output files:
  ## `DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_cophenetic/analysis.collectTrees_iqtree.sh`
  ## output objects:
    1. `/project/berglandlab/alan/cdlo_250K.boot.Rdata` A very large object that has all of the trees, bootstraps, copenhoetic distances
    2. `/project/berglandlab/alan/cdlo_250K.boot.manhattan.Rdata` A small object that has the data to make the genome plot
    3. `/project/berglandlab/alan/cdlo_250K.boot.trees.Rdata` the njo tree for each window. This is only one tree per window, the bootstraps are elsewhere.

# a simple plot
  ## `DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_cophenetic/plot.iqtree.manhattan.R`
  ## `DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_cophenetic/plot.iqtree.tree.R`
