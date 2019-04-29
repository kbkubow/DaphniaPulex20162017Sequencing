First make mpileup of the D8 moderate coverage, no SC a individuals, at the linkage pruned SNPs.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 96:00:00
#SBATCH --mem=64000
#SBATCH -p standard
#SBATCH -A berglandlab

chunk=$1

echo starting chunk ${chunk}

cd /scratch/kbb7sh/Daphnia/NewMapping/mpileups/March2018

echo running samtools

module load samtools 

# Run program

samtools mpileup -q 5 -Q 5 -b /scratch/kbb7sh/Daphnia/NewMapping/mpileups/March2018/March2018noAidsdt \
        -l ${chunk} -o ${chunk}.mpileup
```
Call this script by sbatch mpileup.slurm LDprunedSNPs_20190429.bed.
