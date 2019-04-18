Here is an example of the mapping script used for mapping reads to the HiC, gap closed genome. I included all scaffolds (Daphnia and non-Daphnia) in the reference genome for mapping. But as part of processing the output from BWA, only reads that mapped to Daphnia scaffolds were retained in the resulting bam files.

```
Last login: Thu Apr 18 08:28:29 on ttys001
d-172-25-202-128:~ kbkubow$ ssh -X karen@bergland-lab.bio.virginia.edu
karen@bergland-lab.bio.virginia.edu's password: 
Last failed login: Wed Apr 17 15:22:59 EDT 2019 from 128.143.245.141 on ssh:notty
There was 1 failed login attempt since the last successful login.
Last login: Wed Apr 17 15:21:23 2019 from 137.54.96.6
(base) [karen@bergland-lab ~]$ cd /mnt/spicy_3/Karen/Sp2017/NewMapping
(base) [karen@bergland-lab NewMapping]$ ls
5kbchrassignHiCnew.csv                          LDpruneSNPsBS_A.bed                           RMsnps_20190129.Rdata
Amhetlong_20190221.Rdata                        LDpruneSNPsBS_B.bed                           runingRoH_20190412
computingoutofdpsnps_20190130.Rdata             LDpruneSNPsBS_C.bed                           runingRoH_20190415
consensusAcalls_20181022.Rdata                  LDpruneSNPsBS_D.bed                           runingRoHK10B5_20190415
consensusBcalls_20181022.Rdata                  LDpruneSNPsBS_E.bed                           runingRoHK10B5.R
consensusCcalls_20181022.Rdata                  LDpruneSNPsBS_F.bed                           runingRoHK10B5witherr_20190416
consensusDcalls_20181022.Rdata                  LDpruneSNPsBS_G.bed                           runingRoHK10B5witherr_20190417
consensusEcalls_20181022.Rdata                  LDpruneSNPsBS_H.bed                           runingRoHK10B5witherr.R
consensusFcalls_20181022.Rdata                  LDpruneSNPsBS_I.bed                           runingRoH.Rdata
consensusGcalls_20181022.Rdata                  LDpruneSNPsBS_J.bed                           scsfullnoident_20190131.Rdata
consensusHcalls_20181022.Rdata                  lowhighRDsnps_20190129.Rdata                  sedEavuWB
consensusIcalls_20181022.Rdata                  mABDHIKLMCJGconsensus_20190131.Rdata          sedf0ETzW
consensusJcalls_20181022.Rdata                  MAstousemedrd9                                snpsforbedfinalmap20190220.bed
consensusKcalls_20181022.Rdata                  mat, order = "hclust", method = "color")      subsetindwithcountsctotal_20180821.Rdata
consensusLcalls_20181022.Rdata                  MBstousemedrd9                                superclonesnewmap
consensusMcalls_20181022.Rdata                  MCstousemedrd9                                superclonesnewmap_20190208.Rdata
delly                                           medianandsdreaddepth.Rdata                    test
dpfiltandscfilt_20190131.Rdata                  medianreaddepthsnewsc.Rdata                   testoutfileB.gen
dpfiltsnps_20190130.Rdata                       medianreaddepths.Rdata                        testreal.Rdata
Dpulex_mix10R_B.Rdata                           mito                                          tmp
Dpulex_mix10R.Rdata                             npgdsIBS(genofile, num.thread=2)              totalhamdist_2500_25000_20190204.Rdata
Dpulex_mixK10B5.Rdata                           NsChrRDsnps_20190129.Rdata                    totalhamdist_5000_50000_20190131.Rdata
Dpulex_mixK10B5witherr.Rdata                    numgoodbpsperwindow_20190204.Rdata            totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.gds
finalsetsnpset01_20190131.Rdata                 numgoodbpsperwindowB_20190131.Rdata           totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss_lowDpmiss.gds
finalsnpset_20190131.Rdata                      numgoodbpsperwindowB_20190201.Rdata           totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss_lowDpmiss.seq.gds
fixedwithinpulexids_20190131.Rdata              outfileB.gen                                  totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss_lowDpmiss.vcf
fullPCAbychr_20190131.Rdata                     outfileB.samples                              totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss_lowDpmiss.vcf.idx
fullPCAbychrnoDBarbW_20190131.Rdata             outfile.gen                                   totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.seq.gds
fullPCAbychrnoDBarbWnoD10_20190131.Rdata        outfile.samples                               totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.vcf
fullPCAbychrnoDBarbWnoD10no2016_20190131.Rdata  parseDPGQvcf_20191116                         totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.vcf.idx
gethammingdist_20190131_0131                    parseDPGQvcf.R                                updatesnpstouse_finalssnpssetnorialleleic_20190131AD
gethammingdist_20190131.R                       parseQGvcf_2019116                            updatesnpstouse_finalssnpssetnorialleleic_20190131.Rdata
gethammingdist_25000window_20190131_0204        parseQGvcf.R                                  updatesnpstouse_finalssnpssetnorialleleic_20190131_Samples
gethammingdist_25000window_20190131.R           readdepthrow.ag_medrdperclone_20190130.Rdata  updatesnpstouse_finalssnpssetnorialleleic_20190131.vcf
goodsnpsnotinNsChrRDorRMtable_20190129.Rdata    recalling                                     usebcftoolstoextractAD.sh
goodsnpsnotinNsChrRDtable_20190129.Rdata        Relatedness
(base) [karen@bergland-lab NewMapping]$ pwd
/mnt/spicy_3/Karen/Sp2017/NewMapping
(base) [karen@bergland-lab NewMapping]$ cd /mnt/spicy_3/Karen/scripts/
(base) [karen@bergland-lab scripts]$ ls
Mapping2016update_20180624                 Mapping2017Plate3updatediffpear_D_20181029   MappingHiCGMSpFall2016_20181208  MappingSpFall2016D_090517
Mapping2016update.sh                       Mapping2017Plate3updatediffpear_H_20181029   MappingHiCGMSpFall2016.sh        MappingSpFall2016D84a_090817
Mapping2017P5and6andmakeup_20180508        Mapping2017Plate3updatediffpear_I_20181030   MappingMarch2018SM.sh            MappingSpFall2016D84a_091517
Mapping2017P5and6andmakeup_20180514        Mapping2017Plate3updatediffpear_KL_20181030  MappingMitoP3_20190306           MappingSpFall2016D84a.sh
Mapping2017P5and6andmakeup_20180517        Mapping2017Plate3updatediffpear_M_20181030   MappingMitoP3_20190312           MappingSpFall2016D84Atrial_091517
Mapping2017P5and6andmakeup.sh              Mapping2017Plate3updatediffpear.sh           MappingMitoP3_20190322           MappingSpFall2016D84Atrial.sh
Mapping2017Plate2_20171219                 Mapping2017Plate3update.sh                   MappingMitoP3.sh                 MappingSpFall2016DB_090517
Mapping2017Plate2_20180102                 Mapping2017poolsamp_101317                   MappingMitoP5and6_20190318       MappingSpFall2016DB.sh
Mapping2017Plate2_20180102B                Mapping2017poolsamp.sh                       MappingMitoP5and6.sh             MappingSpFall2016D.sh
Mapping2017Plate2_20180103                 MappingHeli_20180216                         MappingMitoSpFall2016_20190317   MappingSpFall2016E_090417
Mapping2017Plate2_20180105                 MappingHeli.sh                               MappingMitoSpFall2016.sh         MappingSpFall2016E_090517
Mapping2017Plate2_20180105_B               MappingHiCGM_20181203                        Mappingrefgenomereads_20180521   MappingSpFall2016EB_090517
Mapping2017Plate2_20180105C                MappingHiCGM_20181203B                       Mappingrefgenomereads_20180522   MappingSpFall2016EB.sh
Mapping2017Plate2_20180107                 MappingHiCGM_20181203C                       Mappingrefgenomereads.sh         MappingSpFall2016E.sh
Mapping2017Plate2_20180116                 MappingHiCGM_20181204                        mapping.sh                       MappingSpFall2016.sh
Mapping2017Plate2dups_20180116             MappingHiCGM_20181204B                       MappingSp2017D84a_091017         mapReads.sh
Mapping2017Plate2dups.sh                   MappingHiCGM_20181204C                       MappingSp2017D84aB_091017        nohup.out
Mapping2017Plate2.sh                       MappingHiCGM_20181208                        MappingSp2017D84aB_091117        pear_bwa_H2C3CCCXY_0625
Mapping2017Plate3update_20180620           MappingHiCGMP2_20181218                      MappingSp2017D84aB.sh            pear_bwa_H2C3CCCXY_0625B
Mapping2017Plate3updatedHiCnew_20181014    MappingHiCGMP2.sh                            MappingSp2017D84a.sh             pear_bwa_H2C3CCCXY.sh
Mapping2017Plate3updatedHiCnew_20181022    MappingHiCGMP5and6_20181210                  MappingSp2017mito_10162017       pear_bwa_H2CJ5CCXY_0626
Mapping2017Plate3updatedHiCnew_20181023    MappingHiCGMP5and6_20181218                  MappingSp2017mito.sh             pear_bwa_H2CJ5CCXY_0627
Mapping2017Plate3updatedHiCnew_20181023B   MappingHiCGMP5and6.sh                        MappingSpFall2016_080917         pear_bwa_H2CJ5CCXY_0628
Mapping2017Plate3updatedHiCnew_20181024Is  MappingHiCGMPool2012_20191116                MappingSpFall2016B_080917        pear_bwa_H2CJ5CCXY.sh
Mapping2017Plate3updatedHiCnew_20181024KL  MappingHiCGMPool2012.sh                      MappingSpFall2016B.sh            pear_bwa.sh
Mapping2017Plate3updatedHiCnew_20181025    MappingHiCGMPool2017_20190114                MappingSpFall2016C_080917        PriscillaJoinFiles.sh
Mapping2017Plate3updatedHiCnew.sh          MappingHiCGMPool2017.sh                      MappingSpFall2016C.sh
Mapping2017Plate3updatediffpear_20181011   MappingHiCGM.sh                              MappingSpFall2016D_090417
(base) [karen@bergland-lab scripts]$ less MappingHiCGM.sh 

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
