#!/bin/bash

# NOTE: INTERVAL was replaced with a genomic interval of 5000000bp

module load java
module load samtools

chrom=INTERVAL
in_gatk=/Genomics/grid/users/alea/programs/gatk-4.1.4.0
in_genome=/Genomics/ayroleslab2/yushi/ref/hg38_all_chr.fa
in_database=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/$chrom
in_samples=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/all_high_cov.sample_map
out_vcf_gz=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/$chrom.all_high_cov.vcf.gz

# merge GVCFs
$in_gatk/gatk --java-options "-Xmx80g -Xms80g" GenomicsDBImport \
            --sample-name-map $in_samples \
            --genomicsdb-workspace-path $in_database \
            --tmp-dir=/scratch/tmp/ayroles/alea_files \
            -L $chrom \
            --batch-size 40 --reader-threads 16

# call genotypes

echo 'calling joint genotypes...'

echo $chrom > temp.${chrom}.list

$in_gatk/gatk --java-options "-Xmx80g" GenotypeGVCFs \
            -R $in_genome \
            -L temp.${chrom}.list \
            -V gendb://$in_database \
            -O $out_vcf_gz \
            --tmp-dir=/scratch/tmp/ayroles/alea_files --max-alternate-alleles 2
