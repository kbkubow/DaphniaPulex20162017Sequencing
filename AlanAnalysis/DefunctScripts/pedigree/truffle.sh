ijob -c1 -p standard -A berglandlab


### install truffle
  cd /scratch/aob2x/daphnia_hwe_sims/pedigree

  wget https://adimitromanolakis.github.io/truffle-website/downloads/truffle-linux-bin-64-v1.38.tar.gz
  tar xzvf truffle-linux-bin-64-v1.38.tar.gz

### rename chrs
  module load bcftools/1.9
  bcftools \
  annotate \
  --rename-chrs /scratch/aob2x/daphnia_hwe_sims/pedigree/chr_rename.delim \
  -O v \
  -o /scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.vcf \
  /scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.vcf.gz

  zcat /scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.vcf.gz | grep -v "##" | head -n10 | awk '{print NF}'
  cat /scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.vcf | grep -v "##" | head -n10 | awk '{print NF}'


  bcftools \
  view \
  -q 0.01:minor \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.bcf > \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.vcf


### run truffle
  truffle_wd=/scratch/aob2x/daphnia_hwe_sims/pedigree/truffle
  cd $truffle_wd

  ${truffle_wd}/truffle \
  --vcf /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.vcf \
  --cpu 4 \
  --mindist 0 \
  --out /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.vcf.ibd

cp /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.vcf.ibd.ibd \
/nv/vol186/bergland-lab/alan/.


library(data.table)
library(GWASTools)
library(SeqArray)
library(SNPRelate)

ibd <- fread("/mnt/sammas_storage/bergland-lab/alan/12chrs.whatshapp.onePerSC.renameChr.vcf.ibd.ibd")

ibdPlot(k0=ibd$IBD0, k1=ibd$IBD1, relation=ibdAssignRelatedness(k0=ibd$IBD0, k1=ibd$IBD1))


, relation=NULL, color=NULL,
        rel.lwd=2, rel.draw=c("FS", "Deg2", "Deg3"), ...)
