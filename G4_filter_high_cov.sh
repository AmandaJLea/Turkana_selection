#!/bin/sh

module load java
module load samtools
module load bcftools

bcftools concat -a -Oz -o /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.vcf.gz /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/joint_vcfs/*.all_high_cov.vcf.gz

# check for problematic samples
/Genomics/grid/users/alea/programs/plink_1.90 --vcf /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.vcf.gz --snps-only --missing --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/missing

##########
# USE HARD FILTERING FIRST TO REMOVE STUFF WE DEFINITELY WON'T ANALYZE
##########

# region filtering
cpgs=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/CpGs_hg38_v2.bed
mask=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/20160622.allChr.pilot_mask.bed
super=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/genomicSuperDups_hg38.bed
cat $cpgs $super > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp.bed

awk '$4=(FNR FS)' /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp.bed > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v2.bed

awk '$1 !~ /_/' /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v2.bed > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v3.bed

# GATK reccomends filtering for excess hets before VQSR, HWE filtering should do something similar: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
/Genomics/grid/users/alea/programs/plink_1.90 --vcf /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.vcf.gz --snps-only --hwe 0.000001 --make-just-bim --extract range $mask --maf 0.01 --exclude range /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/temp_v3.bed --allow-extra-chr --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs

# turn into bed
awk '{OFS=""; print "chr",$1,"\t",$4-1,"\t",$4}' /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs.bim >  /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs_loc1.bed

awk '{OFS=""; print "chr",$1,"\t",$4,"\t",$4+1}' /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs.bim >  /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs_loc2.bed

tabix -p vcf /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.vcf.gz

# remove low MAF and problem regions, remove 2 problem samples
module load bcftools

regions=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs_loc1.bed
bcftools view -Oz -o /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs1.vcf.gz --samples ^A61,A55 --regions-file $regions /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.vcf.gz

regions=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs_loc2.bed
bcftools view -Oz -o /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs2.vcf.gz --samples ^A61,A55 --regions-file $regions /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.vcf.gz

# filter and sort
bcftools sort -Oz -o /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs1_sort.vcf.gz /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs1.vcf.gz

$in_gatk/gatk IndexFeatureFile -F /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs1_sort.vcf.gz

##########
# USE GATK MIXTURE MODEL FILTERING
##########

echo 'calculating VQSLOD tranches for SNPs...'

in_vcf_gz=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/9May20_high_cov_allCHR.SNPs1_sort.vcf.gz

in_gatk=/Genomics/grid/users/alea/programs/gatk-4.1.4.0
in_genome=/Genomics/ayroleslab2/yushi/ref/hg38_all_chr.fa
in_resourceDIR=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets

out_r_plots=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov_allCHR.SNP1.plots.R
out_snps_recal=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov_allCHR.SNP1.recal
out_snps_tranc=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov_allCHR.SNP1.tranches
out_vcf_gz=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov_allCHR.SNP1.vqsr.vcf.gz

$in_gatk/gatk --java-options "-Xmx60g -Xms60g" VariantRecalibrator \
            -V $in_vcf_gz \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQ -an ReadPosRankSum -an FS -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $in_resourceDIR/hapmap_3.3.hg38.vcf.gz \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 $in_resourceDIR/1000G_omni2.5.hg38.vcf.gz \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 $in_resourceDIR/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
            -O $out_snps_recal \
            --tranches-file $out_snps_tranc \
            -R $in_genome \
            --rscript-file $out_r_plots

echo 'filtering SNPs with VQSLOD...'

$in_gatk/gatk --java-options "-Xmx60g -Xms60g" ApplyVQSR \
            -V $in_vcf_gz \
            --recal-file $out_snps_recal \
            --tranches-file $out_snps_tranc \
            --truth-sensitivity-filter-level 90.0 \
            --create-output-variant-index true \
            -mode SNP \
            -O $out_vcf_gz \
            -R $in_genome

# create one VCF per chromosome
for chr in {1..22}; do /Genomics/grid/users/alea/programs/plink2 --vcf /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov_allCHR.SNP1.vqsr.vcf.gz --var-filter --recode vcf --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov_allCHR.SNP1.${chr}.vcf.gz --chr ${chr}; echo ${chr}; done

