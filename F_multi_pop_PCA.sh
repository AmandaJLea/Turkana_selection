#############
# RUN IN BASH
#############

LDfilt=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_plink
plink=/Genomics/grid/users/alea/programs/plink_1.90
path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
AFR_samples=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/shapeit_files/AFR_1000G_samples.txt

# make bed for 1000G
$plink --vcf $path_resources/shapeit_files/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_w1000G --indep-pairwise 50 20 0.8 --maf 0.01 --keep $AFR_samples

# change rsID
awk '{OFS="\t";print $1,$1":"$4,$3,$4,$5,$6}' /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_w1000G.bim > temp.bim
mv temp.bim /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_w1000G.bim

# first pass merge
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_w1000G --bmerge $LDfilt --out temp1

# filter 1000G
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_w1000G --make-bed --exclude temp1.missnp --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes

# filter Turkana
$plink --bfile $LDfilt --make-bed --exclude temp1.missnp --out temp2

# final merge
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes --make-bed --out /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_w1000G --bmerge temp2 --geno 0.25 --maf 0.01  --indep-pairwise 50 20 0.8

# PCA
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_w1000G --pca 10 tabs --out --geno 0.25 --maf 0.01 --indep-pairwise 50 20 0.8

#############
# RUN IN R
#############

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data")

data1=read.delim('impute1panel.1.merged_LDfilt_w1000G.eigenval.txt',header=F)
data2=read.delim('impute1panel.1.merged_LDfilt_w1000G.eigenvec.txt',header=F)
data3=read.delim('impute1panel.1.merged_LDfilt_w1000G.fam',header=F)
info=read.delim('20130606_g1k.ped.txt')

data2$id<-paste(data2$V1,data2$V2,sep='_')
data3$id<-paste(data3$V1,data3$V2,sep='_')
data4=merge(data2,data3,by='id')
data4=merge(data4,info,by.x='V2.x',by.y='Individual.ID',all.x=T)
# write.table(data4,'impute1panel.1.merged_LDfilt_w1000G.PCA_winfo.txt',row.names=F,sep='\t')

#####
data1=read.delim('impute1panel.1.merged_LDfilt_w1000G.eigenval.txt',header=F)
pve<-c(10.8767,4.5519,2.36929,2.09585,2.07101,2.01789,2.00212,1.962,1.95117,1.93637)

plot(1:10,pve,xlab='PC',ylab='PVE',pch=20,bty='n')
library(RColorBrewer)
cols<-brewer.pal(12,'Paired')
cols2<-c(cols[1:11],'grey',cols[12])

colourCount = length(unique(data4$Population))
getPalette = colorRampPalette(brewer.pal(10, "Paired"))

library(ggplot2)
ggplot(data4, aes(x=V3.x, y=V4.x,col=Population)) +geom_point(alpha=0.65,size=2)+theme_bw(13)+xlab("PC1")+ylab("PC2")+theme_bw(13)+scale_color_manual(values = getPalette(colourCount))
ggplot(data4, aes(x=V5.x, y=V6.x,col=Population)) +geom_point(alpha=0.65,size=2)+theme_bw(13)+xlab("PC3")+ylab("PC4")+theme_bw(13)+scale_color_manual(values = getPalette(colourCount))
ggplot(data4, aes(x=V7, y=V8,col=Population)) +geom_point(alpha=0.65,size=2)+theme_bw(13)+xlab("PC5")+ylab("PC6")+theme_bw(13)+scale_color_manual(values = getPalette(colourCount))


ggplot(data4, aes(x=Population, y=V3.x)) +theme_bw(13)+ylab("PC1")+geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(data4, aes(x=Population, y=V4.x)) +theme_bw(13)+ylab("PC2")+geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

