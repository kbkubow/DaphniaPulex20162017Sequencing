#!/bin/bash

sbatch --array=0-999 \
/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/HWE_simulations.scatter.rivanna.slurm
