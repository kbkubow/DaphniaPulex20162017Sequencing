#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:10:00 # Running time of 4 days
#SBATCH --mem 2G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=1-8 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/coverageDepth/covergeDepth.sh
### sacct -u aob2x -j 20558166
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20451741_1.err

module load samtools gcc/9.2.0  openmpi/3.1.6 python/3.7.7 bedtools/2.29.2

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

start_win=$( expr $start - 1000000 )
stop_win=$( expr $stop + 1000000 )

#cat /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.bed | \
#grep -v "5196681" | grep -v "5201321" > /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.keep.bed

samtools view -b \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam \
${chr}:${start_win}-${stop_win} |
bedtools subtract -A -a - -b /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.keep.bed > \
~/${samp}.small.filter.rg.bam

samtools index ~/${samp}.small.filter.rg.bam








#scp aob2x@rivanna.hpc.virginia.edu:~/*.small.filter.rg.bam* ~/.


#chrLen=$( samtools idxstats /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam | grep ${chr} | cut -f2 )
#chrReads=$( samtools idxstats /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam | grep ${chr} | cut -f3 )

#gzip -c /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa > /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa.gz
#python /scratch/aob2x/dest/DEST/mappingPipeline/scripts/PickleRef.py \
#--ref /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa.gz \
#--output /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.pickle

#samtools mpileup \
#-r ${chr}:${start_win}-${stop_win} \
#-B \
#--excl-flags UNMAP \
#-t DP,AD \
#/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam > ~/foo.txt
#less ~/foo.txt
#q
#module load samtools/0.1.20

#samtools mpileup -AB \
#-r ${chr}:${start_win}-${stop_win} \
#/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam > foo.txt
#cat foo.txt | python /scratch/aob2x/daphnia_hwe_sims/allelecount/allelecount.py


#java -ea -Xmx7g -jar /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/coverageDepth/mpileup2sync.jar \
#--fastq-type sanger \
#--min-qual 0 \
#--input ~/foo.txt \
#--output ~/p1_p2.sync

#less ~/p1_p2.sync

#less -S ~/foo.txt
#less ~/foo.txt | awk -v OFS='\t' '{ if ($4>0 && $5 !~ /[^\^][<>]/ && $5 !~ /\+[0-9]+[ACGTNacgtn]+/ && $5 !~ /-[0-9]+[ACGTNacgtn]+/ && $5 !~ /[^\^]\*/) print $1,$2-1,$2,$3,$4,$5,$6}'

#python /scratch/aob2x/dest/DEST_freeze1/mappingPipeline/scripts/Mpileup2Sync.py \
#--mpileup ~/foo.txt \
#--base-quality-threshold 25 \
#--minIndel 1 \
#--coding 1.8 \
#--ref /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.pickle.ref \
#--output ~/output

#cut -f1,2,4 | sed "s/^/$samp\t${chrLen}\t${chrReads}\t/g" > /scratch/aob2x/daphnia_hwe_sims/rnaseq/coverage/${samp}_${gene}.coveragePos.delim
