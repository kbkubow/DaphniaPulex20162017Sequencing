This is an example of the slurm script I used for making gvcf files for individuals samples using gatk.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 48:00:00
#SBATCH --mem=64000
#SBATCH -p standard
#SBATCH -A berglandlab

chunk=$1

echo starting chunk ${chunk}

cd /scratch/kbb7sh/Daphnia/Sp2017

echo running gatk 

# Run program

java -Xmx4g -jar /scratch/kbb7sh/Daphnia/SpFall2016/GenomeAnalysisTK.jar \
-R /scratch/kbb7sh/Daphnia/NewMapping/totalHiCwithallbestgapclosed.fa \
-T HaplotypeCaller \
-I /scratch/kbb7sh/Daphnia/NewMapping/Plate3HiCGMcloseMap/${chunk}_finalmap.bam \
--emitRefConfidence GVCF \
--variant_index_type LINEAR \
--variant_index_parameter 128000 \
-o /scratch/kbb7sh/Daphnia/NewMapping/newmapgvcf/${chunk}.HiCGM.bam_gvcf
```
