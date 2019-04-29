Let's look at relatedness of the non superclone A March 2018 D8 clones using mapgd. We can use this to test the hypothesis that they are F1 hybrids of A and B, by first looking to see if they are full clones. So just running this analysis on the March2018 individuals for now. Can try adding in A and B later. First make mpileup of the D8 moderate coverage, no SC a individuals, at the linkage pruned SNPs.
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
Call this script by sbatch mpileup.slurm LDprunedSNPs_20190429.bed. Now need to make a header file. 
```
samtools view -H /scratch/kbb7sh/Daphnia/SingleMomsMarch2018/March2018SMMap/March20_2018_D8_1.sort.mdup.bam > seq1.header
```
Next, generate a pro file.
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

echo running mapgd

module load singularity 

# Run program

singularity run ./../mapgd-srf_latest.sif proview -i ${chunk}.mpileup -H seq1.header -s > ${chunk}.pro;
```
Call this script by sbatch pro.slurm LDprunedSNPs_20190429.bed. Next make a map file.
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

echo running mapgd

module load singularity 

# Run program

singularity run ./../mapgd-srf_latest.sif allele -i ${chunk}.pro -o ${chunk}.map;
```
Call this script by sbatch map.slurm LDprunedSNPs_20190429.bed.
Now filter the map file.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 72:00:00
#SBATCH --mem=64000
#SBATCH -p standard
#SBATCH -A berglandlab

chunk=$1

echo starting chunk ${chunk}

cd /scratch/kbb7sh/Daphnia/NewMapping/mpileups/March2018

echo running mapgd

module load singularity 

# Run program

singularity run ./../mapgd-srf_latest.sif filter -i ${chunk}.map.map -p 22 -E 0.01 -c 5000 -C 8000 > ${chunk}.filter.map
```
Call this script by sbatch filter.slurm LDprunedSNPs_20190429.bed.
Then create a genotype file.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 72:00:00
#SBATCH --mem=64000
#SBATCH -p standard
#SBATCH -A berglandlab

chunk=$1

echo starting chunk ${chunk}

cd /scratch/kbb7sh/Daphnia/NewMapping/mpileups/March2018

echo running mapgd

module load singularity 

# Run program

singularity run ./../mapgd-srf_latest.sif genotype --map ${chunk}.filter.map --pro ${chunk}.pro -o ${chunk};
```
Call this script by sbatch geno.slurm LDprunedSNPs_20190429.bed.
Now try running relatedness.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 72:00:00
#SBATCH --mem=250000
#SBATCH -p largemem
#SBATCH -A berglandlab

chunk=$1

echo starting chunk ${chunk}

cd /scratch/kbb7sh/Daphnia/NewMapping/mpileups/March2018

echo running mapgd

module load singularity 

# Run program

singularity run ./../mapgd-srf_latest.sif relatedness -l -i ${chunk}.gcf -o ${chunk}relate
```
Call this script by sbatch relate.slurm LDprunedSNPs_20190429.bed.







