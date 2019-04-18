Filtering my vcf.
This is a redo of an earlier filtering, where I removed all snps within 25 basepairs of indels and set low read depth genotype calls to missing in addition to low genotype quality genotypes. I am now changing this to removing snps within 10 basepairs of indels, and only setting low quality genotype scores to missing.
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

cd /scratch/kbb7sh/Daphnia/NewMapping/finalvcf

bcftools filter --SnpGap 10 totalnewmapB.vcf -o totalnewmapBfiltsnps10bpindels.vcf
```
Line count of original file was 7,885,686. Line count of ouput file was 5,694,221. Assuming this change in line count accurately reflects the number of snps removed, then 2,191,465 snps were removed.
Next remove all indels, so the file has only SNPs.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 4:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/NewMapping/finalvcf

echo running gatk 

java -Xmx4g -jar /scratch/kbb7sh/Daphnia/SpFall2016/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R /scratch/kbb7sh/Daphnia/NewMapping/totalHiCwithallbestgapclosed.fa \
	-V /scratch/kbb7sh/Daphnia/NewMapping/finalvcf/totalnewmapBfiltsnps10bpindels.vcf \
    -selectType SNP \
    -o totalnewmapBfiltsnps10bpindels_snps.vcf 
```
Line count of file after indels are removed is 3,243,021. Assuming this change in line count accurately reflects the number of indels removed, then 2,451,1200 indels were removed.
Next hard filter the SNPs based on GATK recommendations for organisms with no reference SNP panel.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 4:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/NewMapping/finalvcf

echo running gatk 

# Run program

java -Xmx4g -jar /scratch/kbb7sh/Daphnia/SpFall2016/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R /scratch/kbb7sh/Daphnia/NewMapping/totalHiCwithallbestgapclosed.fa \
	-V /scratch/kbb7sh/Daphnia/NewMapping/finalvcf/totalnewmapBfiltsnps10bpindels_snps.vcf   \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "my_snp_filter" \
    -o totalnewmapBfiltsnps10bpindels_snps_filter.vcf
---
---
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/NewMapping/finalvcf

echo running gatk 

# Run program


java -Xmx4g -jar /scratch/kbb7sh/Daphnia/SpFall2016/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R /scratch/kbb7sh/Daphnia/NewMapping/totalHiCwithallbestgapclosed.fa \
	-V /scratch/kbb7sh/Daphnia/NewMapping/finalvcf/totalnewmapBfiltsnps10bpindels_snps_filter.vcf   \
	-ef \
    -o totalnewmapBfiltsnps10bpindels_snps_filter_pass.vcf
```
Filtering out SNPs results in a file with a line count of 2,862,375. Assuming the change in line count reflects the change in SNPs, then 380,646 SNPs were dropped based on filtering.
Set low GQ (scores less than 10) to zero.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/NewMapping/finalvcf

module load bcftools

# Run program

bcftools filter -e "FORMAT/GQ<10" -S "." totalnewmapBfiltsnps10bpindels_snps_filter_pass.vcf | bcftools view -O v -o totalnewmapBfiltsnps10bpindels_snps_filter_pass_lowGQmiss.vcf
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

cd /scratch/kbb7sh/Daphnia/NewMapping/finalvcf

 java -Xmx4g -jar /scratch/kbb7sh/Daphnia/SpFall2016/GenomeAnalysisTK.jar \
   -T ValidateVariants \
   -R /scratch/kbb7sh/Daphnia/NewMapping/totalHiCwithallbestgapclosed.fa \
   -V totalnewmapBfiltsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
   --validationTypeToExclude ALL
```
The final filtered VCF has 2,853,063 SNP variants. Next step will be to import the VCF into R and do further filtering.
