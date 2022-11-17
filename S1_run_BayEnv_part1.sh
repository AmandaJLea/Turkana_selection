#!/bin/bash

# NOTE: CHROMNUM was replaced with the chromosome number

#######
# high and low/med coverage WGS data, LD filtered - add LWK and YRI data
#######

chr=CHROMNUM
LDfilt=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_plink
plink=/Genomics/grid/users/alea/programs/plink_1.90
path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
AFR_samples=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/ALL_YRI_LWK_samples.txt

# make bed for 1000G
$plink --vcf $path_resources/shapeit_files/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G --indep-pairwise 50 20 0.8 --maf 0.01 --keep $AFR_samples

# change rsID
awk '{OFS="\t";print $1,$1":"$4,$3,$4,$5,$6}' /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G.bim > temp_${chr}.bim
mv temp_${chr}.bim /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G.bim

# first pass merge
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G --bmerge $LDfilt --out temp1_${chr}

# filter 1000G
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G --make-bed --exclude temp1_${chr}.missnp --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes

# filter Turkana
$plink --bfile $LDfilt --make-bed --exclude temp1_${chr}.missnp --out temp2_${chr}

# final merge
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes --make-bed --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G --bmerge temp2_${chr} --geno 0.25 --maf 0.01  --indep-pairwise 50 20 0.8

#######
# calculate allele frequency by group
#######

# cat /Genomics/ayroleslab2/alea2/16May20_samples_to_keep.txt /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/ALL_YRI_LWK_samples.txt > /Genomics/ayroleslab2/alea2/16May20_samples_to_keep_wYRILWK.txt 

plink=/Genomics/grid/users/alea/programs/plink_1.90
groups=/Genomics/ayroleslab2/alea2/16May20_groups_wYRILWK.txt
keep=/Genomics/ayroleslab2/alea2/16May20_samples_to_keep_wYRILWK.txt
groups=/Genomics/ayroleslab2/alea2/16May20_groups_wYRILWK.txt

$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G --extract /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G.prune.in --make-bed --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt2_w1000G --keep $keep --geno 0.5 --maf 0.05 --freq --within $groups

groups=/Genomics/ayroleslab2/alea2/16May20_groups_wjust_TURKYRILWK.txt
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G --extract /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt_w1000G.prune.in --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.${chr}.merged_LDfilt3_w1000G --keep $keep --geno 0.25 --maf 0.01 --fst --within $groups

