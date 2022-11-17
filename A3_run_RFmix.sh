############
# RFMix info
############

A phased VCF/BCF file containing "query" haplotypes which are to be analyzed.
A phased VCF/BCF file containing reference haplotypes (in any order)
A reference sample map file matching reference samples to their respective reference populations
A genetic map file
All files mapped or referenced to the same genome assembly
All diploid genotypes must be phase resolved
In addition to the above, you must at minimum also specify a basename (prefix) for output files, and the chromosome to analyze from the VCF/BCF inputs, even if they contain only one chromosome
If BCF files are used, bcftools must be installed and available in the PATH environment setting
VCF/BCF files may be gzip compressed, and should be indexed using bcftools

RFMIX upon completion will output two main files of interest: the most likely assignment of subpopulations per CRF point (<output basename>.msp.tsv), and the marginal probabilities of each subpopulation being the ancestral population of the corresponding CRF point (<output basename>.fb.tsv). 

#############
# RUN IN BASH
# NOTE: CHROM is replaced with the chromosome name
#############

#!/bin/bash

f=CHROM

path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
phased_1KG_haps=$path_resources/shapeit_files/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

# subset 1000 Genomes to LWK and CEU
module load bcftools
bcftools view -O z -o $path_resources/shapeit_files/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.YRI_LWK_CEU.vcf.gz --samples-file YRI_LWK_CEU.txt $phased_1KG_haps --force-samples

# convert to vcf
phased_haps=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov.SNP1.hg19_chr${f}.phased
phased_haps2=high_cov.SNP1.hg19_chr${f}.phased.vcf

/Genomics/grid/users/alea/programs/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit -convert \
        --input-haps $phased_haps \
        --output-vcf $phased_haps2 --thread=4

# run RFmix
# map_path=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_chr${f}_combined_b37.txt
# awk -v chr2=${chr} '{OFS="\t"; print chr2,$1,$3}' $map_path | tail -n +2 > /Genomics/ayroleslab2/alea/tsimane/int_data_files/temp.genetic_map_chr${chr}.txt

map_path=/Genomics/ayroleslab2/alea/tsimane/int_data_files/temp.genetic_map_chr${f}.txt

~/programs/rfmix-master/rfmix -f high_cov.SNP1.hg19_chr${f}.phased.vcf -r $path_resources/shapeit_files/ALL.chr${f}.phase3_shapeit2_mvncall_integrated_v5a.20130502.YRI_LWK_CEU.vcf.gz -m rfmix_info.txt -g $map_path -o RFmix_chr${f} --chromosome=${f} -G 10 

#############
# RUN IN R
#############

# concat RFmix proportions
library(data.table)

#  the most likely assignment of subpopulations per CRF point 
data1_msp=fread('RFmix_chr1.msp.tsv',skip=1)

for (f in 2:22) {
data1_part<-fread(paste('RFmix_chr',f,'.msp.tsv',sep=''),skip=1)
data1_msp=rbind(data1_msp,data1_part) }
data1_msp$count_0<-apply(data1_msp[,c(7:222),with=F],1,function(x) length(which(x==0)))
data1_msp$count_1<-apply(data1_msp[,c(7:222),with=F],1,function(x) length(which(x==1)))
data1_msp$count_2<-apply(data1_msp[,c(7:222),with=F],1,function(x) length(which(x==2)))
write.table(data1_msp[,c(1:6,223:225),with=F],'20Aug22_RFmix_summary.txt',row.names=F,sep='\t')

 
# admixture proportions per chromosome
data1_Q=fread('RFmix_chr1.rfmix.Q',skip=1)

for (f in 2:22) {
data1_part<-fread(paste('RFmix_chr',f,'.rfmix.Q',sep=''),skip=1)
data1_Q=rbind(data1_Q,data1_part) }

write.table(data1_Q,'20Aug22_RFmix_summary2.txt',row.names=F,sep='\t')

# compare ADMIXTURE vs RFmix proportions
setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures")
data1=read.delim('20Aug22_RFmix_vs_ADMIXTURE.txt')

par(mfrow=c(1,2))
plot(data1$CEU,data1$V2,xlab='European ancestry proportion (RFMix)',ylab='European ancestry proportion (ADMIXTURE)',xlim=c(0,0.2),ylim=c(0,0.2),bty='n')
cor.test(data1$CEU,data1$V2)

plot(density(data1$CEU),main='',bty='n',lwd=2,xlab='European ancestry proportion')
lines(density(data1$V2,na.rm=T),lty=2,lwd=2)

