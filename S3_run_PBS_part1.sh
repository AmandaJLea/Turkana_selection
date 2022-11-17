#!/bin/bash

# NOTE: CHROMNUM was replaced with the chromosome number

#######
# high and low/med coverage WGS data, LD filtered - add LWK and YRI data
#######

LDfilt=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_plink
plink=/Genomics/grid/users/alea/programs/plink_1.90
path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
AFR_samples=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/shapeit_files/YRI_LWK_samples.txt

# make bed for 1000G
$plink --vcf $path_resources/shapeit_files/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out impute1panel.${chr}.merged_LDfilt_w1000G --indep-pairwise 50 20 0.8 --maf 0.01 --keep $AFR_samples

# change rsID
awk '{OFS="\t";print $1,$1":"$4,$3,$4,$5,$6}' impute1panel.${chr}.merged_LDfilt_w1000G.bim > temp_${chr}.bim
mv temp_${chr}.bim impute1panel.${chr}.merged_LDfilt_w1000G.bim

# first pass merge
$plink --bfile impute1panel.${chr}.merged_LDfilt_w1000G --bmerge $LDfilt --out temp1_${chr}

# filter 1000G
$plink --bfile impute1panel.${chr}.merged_LDfilt_w1000G --make-bed --exclude temp1_${chr}.missnp --out ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes

# filter Turkana
$plink --bfile $LDfilt --make-bed --exclude temp1_${chr}.missnp --out temp2_${chr}

# final merge
$plink --bfile ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes --make-bed --out impute1panel.${chr}.merged_LDfilt_w1000G --bmerge temp2_${chr} --geno 0.25 --maf 0.01  --indep-pairwise 50 20 0.8

#######
# calculate Fst for each pop comparison
#######

grep 'Turk' LWK_YRI_Turk.groups.txt > YRI_Turk.txt
grep 'YRI' LWK_YRI_Turk.groups.txt >> YRI_Turk.txt

grep 'Turk' LWK_YRI_Turk.groups.txt > LWK_Turk.txt
grep 'LWK' LWK_YRI_Turk.groups.txt >> LWK_Turk.txt

grep 'YRI' LWK_YRI_Turk.groups.txt > LWK_YRI.txt
grep 'LWK' LWK_YRI_Turk.groups.txt >> LWK_YRI.txt


$plink --bfile impute1panel.${chr}.merged_LDfilt_w1000G --out YRI_Turk_chr${chr}_FST --extract impute1panel.${chr}.merged_LDfilt_w1000G.prune.in --within YRI_Turk.txt --fst

$plink --bfile impute1panel.${chr}.merged_LDfilt_w1000G --out LWK_Turk_chr${chr}_FST --extract impute1panel.${chr}.merged_LDfilt_w1000G.prune.in --within LWK_Turk.txt --fst

$plink --bfile impute1panel.${chr}.merged_LDfilt_w1000G --out LWK_YRI_chr${chr}_FST --extract impute1panel.${chr}.merged_LDfilt_w1000G.prune.in --within LWK_YRI.txt --fst

$plink --bfile impute1panel.${chr}.merged_LDfilt_w1000G --out LWK_YRI_Turk_chr${chr} --extract impute1panel.${chr}.merged_LDfilt_w1000G.prune.in --within LWK_YRI_Turk.groups.txt --freq

rm ALL.chr${chr}*
rm temp1_${chr}*
rm temp2_${chr}*

awk '{OFS="\t"; print $2,$3,$4,$5,$6,$8}' LWK_YRI_Turk_chr${chr}.frq.strat  > LWK_YRI_Turk_chr${chr}.frq_v2.strat

