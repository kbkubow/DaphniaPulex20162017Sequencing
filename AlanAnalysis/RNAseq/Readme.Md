### 1. Basic quality control
  `DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/rna_seq_fastqc.sh`

  ```
  module load gcc/9.2.0  openmpi/3.1.6 python/3.7.7
  /home/aob2x/.local/bin/multiqc -f .
  tar czv multiqc_data/ multiqc_report.html > ~/daphnia_rnaseq.multiqc.tar.gz
  ```
  scp aob2x@rivanna.hpc.virginia.edu:~/daphnia_rnaseq.multiqc.tar.gz ~/.

### 2. Adapter removal & paired end merging
