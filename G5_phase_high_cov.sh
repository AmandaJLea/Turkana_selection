#!/bin/sh

# NOTE: NUMBER was replaced with the chromosome number

module load java
module load samtools
module load bcftools

chr=NUMBER

# liftover to hg19
in_picard=/Genomics/grid/users/alea/programs/picard-tools-2.21.6/picard.jar
in_chain=/Genomics/grid/users/alea/programs/hg38ToHg19_nochr.over.chain
in_genome=/Genomics/ayroleslab2/yushi/ref/hg37_10kb/hg19_v0_Homo_sapiens_assembly19.fasta

# define output
out_reject=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.reject_chr${chr}.vcf
out_lifted_vcf=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.vcf

grep -v '*' /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov_allCHR.SNP1.${chr}.vcf.gz.vcf > temp_${chr}.vcf
in_vcf=temp_${chr}.vcf

java -jar $in_picard LiftoverVcf \
   I=$in_vcf \
   OUTPUT=$out_lifted_vcf \
   CHAIN=$in_chain \
   REJECT=$out_reject \
   R=$in_genome WARN_ON_MISSING_CONTIG=true

java -jar $in_picard SortVcf \
      I=$out_lifted_vcf O=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.sort.vcf

/Genomics/grid/users/alea/programs/plink2 --recode vcf bgz --vcf /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.sort.vcf --snps-only --max-alleles 2 --min-alleles 2 --maf 0.01 --keep-allele-order --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.sort --geno 0.05 --chr ${chr}

# phase 
map=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.txt
phased_haps=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.phased.haps
phased_sample=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.phased.sample

/Genomics/grid/users/alea/programs/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit --input-vcf /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.sort.vcf.gz \
        --input-map $map \
        --output-max $phased_haps $phased_sample 

# prep for impute2
phased_haps=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.phased
phased_haps2=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.phased_v2.haps
phased_leg=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.phased_v2.leg
phased_sam=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.phased_v2.sam

/Genomics/grid/users/alea/programs/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit -convert \
        --input-haps $phased_haps \
        --output-ref $phased_haps2 $phased_leg $phased_sam --thread=4
