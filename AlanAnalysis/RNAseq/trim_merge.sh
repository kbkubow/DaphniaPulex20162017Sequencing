#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-01:00:00 # Running time of 4 days
#SBATCH --mem 5G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/rna_seq_fastq.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/rna_seq_fastq.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/rna_seq_fastqc.sh
### sacct -u aob2x -j 20373409
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/rna_seq_fastq.20373409_

#SLURM_ARRAY_TASK_ID=1
samp=$( sed "${SLURM_ARRAY_TASK_ID};d" DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/samples)


~/adapterremoval-2.3.1/build/AdapterRemoval \
