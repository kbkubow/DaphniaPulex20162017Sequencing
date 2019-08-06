#!/bin/bash

sbatch --array=0-501 \
/scratch/aob2x/daphnia_hwe_sims/HW_simulations.scatter.rivanna.slurm
