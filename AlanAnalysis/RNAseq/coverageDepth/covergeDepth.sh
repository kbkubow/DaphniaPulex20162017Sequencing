#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:30:00 # Running time of 4 days
#SBATCH --mem 10G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=1-8 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/coverageDepth/covergeDepth.sh
### sacct -u aob2x -j 20498639
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20451741_1.err

module load samtools

#SLURM_ARRAY_TASK_ID=8
wd=/scratch/aob2x/daphnia_hwe_sims/

samp=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/samples )
echo $samp

if [ ! -f /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam.bai ]; then
  samtools index /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam
fi

gene=Daphnia00787

chr=$( grep   ${gene} /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff | grep -E "mRNA" | cut -f1 )
start=$( grep ${gene} /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff | grep -E "mRNA" | cut -f4 )
stop=$( grep  ${gene} /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff | grep -E "mRNA" | cut -f5 )

start_win=$( expr $start - 10000 )
stop_win=$( expr $start + 10000 )

chrLen=$( samtools idxstats /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam | grep ${chr} | cut -f1 )
chrReads=$( samtools idxstats /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam | grep ${chr} | cut -f2 )

samtools mpileup \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam \
-r ${chr}:${start_win}-${stop_win} | cut -f1,2,4 | sed "s/^/$samp\t${chrLen}\t${chrReads}\t/g" > /scratch/aob2x/daphnia_hwe_sims/rnaseq/coverage/${samp}_${gene}.coveragePos.delim
