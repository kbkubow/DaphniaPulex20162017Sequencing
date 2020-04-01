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



### run truffle
  truffle_wd=/scratch/aob2x/daphnia_hwe_sims/pedigree/truffle
  cd $truffle_wd

  ${truffle_wd}/truffle \
  --vcf /scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.vcf  \
  --cpu 4 \
  --maf 0.01 \
  --missing 0.05 \
  --mindist 0 \
  --out /scratch/aob2x/daphnia_hwe_sims/pedigree/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd
