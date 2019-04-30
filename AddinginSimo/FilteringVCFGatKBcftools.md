Filtering my vcf that includes individual genome sequences from 2016,2017, and March 2018, now also adding in 2 Simo. Also included 2 libraries made from single males of known genotypes.
I am removing snps within 10 basepairs of indels, hard filtering according to gatk's recommendations, and setting low quality genotype scores to missing.
First remove all snps within 10 basepairs of indels.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 1:00:00
#SBATCH --mem=60000
#SBATCH -p standard
#SBATCH -A berglandlab

module load bcftools

cd /scratch/kbb7sh/Daphnia/SingleMomsMarch2018

bcftools filter --SnpGap 10 totalnewmapwMarch2018_D.vcf -o totalnewmapwMarch2018_Dfiltsnps10bpindels.vcf
```
Line count of original file was 9,745,081. Line count of ouput file was 6,757,744. Assuming this change in line count accurately reflects the number of snps removed, then 2,987,337 snps were removed.
Next remove all indels, so the file has only SNPs.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 4:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/SingleMomsMarch2018

echo running gatk 

java -Xmx4g -jar /scratch/kbb7sh/Daphnia/SpFall2016/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R /scratch/kbb7sh/Daphnia/NewMapping/totalHiCwithallbestgapclosed.fa \
	-V /scratch/kbb7sh/Daphnia/SingleMomsMarch2018/totalnewmapwMarch2018_Dfiltsnps10bpindels.vcf \
    -selectType SNP \
    -o totalnewmapwMarch2018_Dfiltsnps10bpindels_snps.vcf 
```
Line count of file after indels are removed is 3,694,923. Assuming this change in line count accurately reflects the number of indels removed, then 2,940,692 indels were removed.
Next hard filter the SNPs based on GATK recommendations for organisms with no reference SNP panel.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 4:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/SingleMomsMarch2018

echo running gatk 

# Run program

java -Xmx4g -jar /scratch/kbb7sh/Daphnia/SpFall2016/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R /scratch/kbb7sh/Daphnia/NewMapping/totalHiCwithallbestgapclosed.fa \
	-V /scratch/kbb7sh/Daphnia/SingleMomsMarch2018/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps.vcf   \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "my_snp_filter" \
    -o totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter.vcf
---
---
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/SingleMomsMarch2018

echo running gatk 

# Run program


java -Xmx4g -jar /scratch/kbb7sh/Daphnia/SpFall2016/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R /scratch/kbb7sh/Daphnia/NewMapping/totalHiCwithallbestgapclosed.fa \
	-V /scratch/kbb7sh/Daphnia/SingleMomsMarch2018/totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter.vcf   \
	-ef \
    -o totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass.vcf
```
Filtering out SNPs results in a file with a line count of 2,998,158. Assuming the change in line count reflects the change in SNPs, then 696,765 SNPs were dropped based on filtering.
Set low GQ (scores less than 10) to zero.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/SingleMomsMarch2018

module load bcftools

# Run program

bcftools filter -e "FORMAT/GQ<10" -S "." totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass.vcf | bcftools view -O v -o totalnewmapwMarch2018_Afiltsnps10bpindels_snps_filter_pass_lowGQmiss.vcf
```
Use gatk check variants to make an index.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/SingleMomsMarch2018

 java -Xmx4g -jar /scratch/kbb7sh/Daphnia/SpFall2016/GenomeAnalysisTK.jar \
   -T ValidateVariants \
   -R /scratch/kbb7sh/Daphnia/NewMapping/totalHiCwithallbestgapclosed.fa \
   -V totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
   --validationTypeToExclude ALL
```
The final filtered VCF has 2,988,846 SNP variants. Next step will be to import the VCF into R and do further filtering.