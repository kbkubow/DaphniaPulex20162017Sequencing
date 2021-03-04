java -Xmx1G -jar ~/QoRTs.jar \
makeFlatGff \
--DEXSeqFmt \
/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf \
/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.flatgff


java -Xmx1G -jar ~/QoRTs.jar \
mergeNovelSplices \
--minCount 6 \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/ \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/decoder.txt \
/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/
