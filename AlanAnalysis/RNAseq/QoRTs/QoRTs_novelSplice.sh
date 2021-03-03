java -Xmx1G -jar ~/QoRTs.jar \
makeFlatGff \
--DEXSeqFmt \
/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf \
/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.flatgff


java -Xmx1G -jar
softwareRelease/QoRTs.jar \
mergeNovelSplices \
--minCount 6 \
outputData/countTables/ \
outputData/sizeFactors.GEO.txt \
inputData/annoFiles/anno.gtf.gz \
outputData/countTables/
