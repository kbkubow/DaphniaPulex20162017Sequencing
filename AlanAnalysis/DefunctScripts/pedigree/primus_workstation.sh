
cp /mnt/sammas_storage/bergland-lab/alan/PRIMUS_v1.9.0.tgz /mnt/internal_1/primus/.
cp /mnt/sammas_storage/bergland-lab/alan/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim /mnt/internal_1/primus/.

wd=/mnt/internal_1/primus/PRIMUS_v1.9.0/bin/

#rm -fr /mnt/internal_1/primus/primusOut/*
cd /mnt/internal_1/primus/primusOut/

nohup ${wd}/run_PRIMUS.pl \
-i FILE=/mnt/internal_1/primus/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim \
IBD0=7 IBD1=8 IBD2=9 PI_HAT=10 \
--max_gen_gap 2 \
-o /mnt/internal_1/primus/primusOut/ \
-v 3 &


scp bergland@bergland-lab.bio.virginia.edu:/mnt/internal_1/primus/primusOut/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim_network2.dot .


scp bergland@bergland-lab.bio.virginia.edu:/mnt/internal_1/primus/primusOut/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim_network1/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim_network1_5.ps ~/.
