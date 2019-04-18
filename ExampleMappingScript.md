Here is an example of the mapping script used for mapping reads to the HiC, gap closed genome. I included all scaffolds (Daphnia and non-Daphnia) in the reference genome for mapping. But as part of processing the output from BWA, only reads that mapped to Daphnia scaffolds were retained in the resulting bam files.

```
#!/bin/bash

# let us stop the whole script on CTRL+C
        int_handler()
        {
                echo "Interrupted."
                # Kill the parent process of the script.
                kill $PPID
                exit 1
        }
        trap 'int_handler' INT

### functions

        
        mapReads () {

                ## define some parameters
                sample=${1}
                sampleB=($(echo ${sample} | cut -d"_" -f1-4))
                sampID=($(echo ${sample} | cut -d"_" -f3))
                pond=($(echo ${sample} | cut -d"_" -f5-8))
                        
                threads=10
                echo $sample
                echo $sampleB
                inputDir="/mnt/spicy_3/Karen/Sp2017/Plate3"
                interDir="/mnt/spicy_3/Karen/Sp2017/Plate3HiCGMcloseInt"
                outputDir="/mnt/spicy_3/Karen/Sp2017/Plate3HiCGMcloseMap"
                        
                flowcell=($( ls ${inputDir}/*.gz | awk '{split($0,a,"/"); print a[7]}' | awk '{split($0,a,"_"); print a[1]}' | sort | uniq ))
                echo ${!flowcell[@]}            
                                
                        ## map reads
                        for cell in "${flowcell[@]}"; do
                                
                                echo "${cell}"
                                lanes=($( ls ${inputDir}/${cell}*${sampleB}* | grep -oE 's[1-8]{1,}' | sort | uniq ))
                                        
                                        for lane in "${lanes[@]}"; do
                                                
                                               echo "${lane}"
                                                ### trim out nextera and index seq
                                                java -jar /mnt/spicy_3/Karen/Trimmomatic-0.36/trimmomatic-0.36.jar PE ${inputDir}/${cell}_${lane}_1_${sampleB}.fastq.gz ${inputDir}/${cell}_${lane}_2_${sampleB}.fastq.gz ${interDir}/${cell}_${lane}_1_${sampID}.P_trimm.fastq ${interDir}/${cell}_${lane}_1_${sampID}.U_trimm.fastq ${interDir}/${cell}_${lane}_2_${sampID}.P_trimm.fastq ${interDir}/${cell}_${lane}_2_${sampID}.U_trimm.fastq ILLUMINACLIP:/mnt/spicy_3/Karen/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:8:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
                
                                                 ### next, map to reference genome
                                                        bwa mem -t ${threads} -K 100000000 -Y \
                                                                -R "@RG\tID:${sampID};${cell};${lane}\tSM:${pond}\tPL:illumina\tPU:${sampID};${cell};${lane}" \
                                                                /mnt/spicy_3/Karen/RefGenome/Dovetail/HiCnew/totalHiCwithallbestgapclosed.fa \
                                                                 ${interDir}/${cell}_${lane}_1_${sampID}.P_trimm.fastq \
                                                                 ${interDir}/${cell}_${lane}_2_${sampID}.P_trimm.fastq | \
                                                        samtools view -L /mnt/spicy_3/Karen/RefGenome/Dovetail/HiCnew/D84Agoodscaffstouse.bed -Suh - | \
                                                        samtools sort -@ ${threads} -o ${interDir}/${cell}_${lane}_${pond}.sort.bam
                                                        samtools index ${interDir}/${cell}_${lane}_${pond}.sort.bam
                                                        java -jar /usr/local/bin/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=${interDir}/${cell}_${lane}_${pond}.sort.bam OUTPUT=${outputDir}/${cell}_${lane}_${pond}.sort.mdup.bam METRICS_FILE=${interDir}/${cell}_${lane}_${pond}.sort.mdups.metrics CREATE_INDEX=true
                                                        
                                                ### next, clean up files in INT directory
                                                rm ${interDir}/*.fastq
                                                rm ${interDir}/*.bam*

                                        done

                        done

                                                ### next, merge bam files to single bam file
                                                samtools merge ${outputDir}/${pond}_finalmap.bam ${outputDir}/*${pond}.sort.mdup.bam
                                                samtools index ${outputDir}/${pond}_finalmap.bam        

                                                ### finally, remove individual bam files
                                                rm ${outputDir}/${cell}*
                                                 
       
        }
        export -f mapReads

### mapping
        parallel --gnu -j1 -a /mnt/spicy_3/Karen/Sp2017/Plate3/inputsnoblanksB mapReads 
```
