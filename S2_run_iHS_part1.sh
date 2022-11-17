#!/bin/bash

# NOTE: chr and CHROMNUM were replaced with the chromosome number

######
# WGS data, high coverage, phased haplotypes - convert from haplotype to VCF format
######

phased_haps=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${chr}.phased
phased_haps_vcf=/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/rehh_files/high_cov.SNP1.hg19_chr${chr}.phased.vcf

/Genomics/grid/users/alea/programs/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit -convert \
        --input-haps $phased_haps \
        --output-vcf $phased_haps_vcf --thread=4

#######
# add ancestral allele info
#######

# wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2

cd /Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/rehh_files

export PERL5LIB=$PERL5LIB:~/programs/vcftools/src/perl/
module load vcftools
module load samtools
phased_haps_vcf=/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/rehh_files/high_cov.SNP1.hg19_chrCHROMNUM.phased.vcf

sed 's,^>.*,>CHROMNUM,' /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/human_ancestor_GRCh37_e59/human_ancestor_CHROMNUM.fa | bgzip > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/human_ancestor_GRCh37_e59/human_ancestor_CHROMNUM.fa.gz
samtools faidx /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/human_ancestor_GRCh37_e59/human_ancestor_CHROMNUM.fa.gz

cat $phased_haps_vcf | fill-aa -a /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/human_ancestor_GRCh37_e59/human_ancestor_CHROMNUM.fa.gz | bgzip -c > high_cov.hg19_chrCHROMNUM.phased_withAA.vcf.gz







