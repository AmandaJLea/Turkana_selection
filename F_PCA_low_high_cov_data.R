###################
# RUN IN BASH
# PCA in Plink
###################

# LD filtered, merged set of high and low coverage samples
LDfilt=/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_plink
plink=/Genomics/grid/users/alea/programs/plink_1.90
path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
AFR_samples=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/shapeit_files/AFR_1000G_samples.txt

$plink --bfile $LDfilt --out $LDfilt --pca 10 tabs

###################
# RUN IN R
# Plot PCA
###################

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data")

data1=read.delim('impute1panel.1.merged_LDfilt_plink.eigenval',header=F)
data2=read.delim('impute1panel.1.merged_LDfilt_plink.eigenvec',header=F)
data3=read.delim('impute1panel.1.merged_LDfilt_plink.fam_winfo.txt',header=F)

data2$id<-paste(data2$V1,data2$V2,sep='_')
data3$id<-paste(data3$V1,data3$V2,sep='_')
data4=merge(data2,data3,by='id')
data4$info<-paste(data4$V3.y,data4$V4.y,sep='_')

plot(1:10,data1$V1,xlab='PC',ylab='PVE',pch=20,bty='n',cex.axis=1.2,cex.lab=1.2)

library(ggplot2)
ggplot(data4, aes(x=V3.x, y=V4.x, col=info,shape=V5.y)) +geom_point(alpha=0.65,size=2)+theme_bw(15)+xlab("PC 1")+ylab("PC 2")+ scale_color_brewer(palette="Paired")

ggplot(data4, aes(x=info, y=V3.x,col=info)) +theme_bw(13)+ylab("PC 1")+geom_violin(width=2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_brewer(palette="Paired")+geom_boxplot(width=0.2)

ggplot(data4, aes(x=info, y=V4.x,col=info)) +theme_bw(13)+ylab("PC 2")+geom_violin(width=2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_brewer(palette="Paired")+geom_boxplot(width=0.2)+ylim(-0.1,0.1)

# ASU vs Princeton
tmp<-subset(data4,V5.y=='high')
pval1<-c()
for (i in 1:10){
pval1<-c(pval1,wilcox.test(tmp[,3+i]~tmp$V4.y)$p.value)}

# low vs high
tmp<-subset(data4,V3.y=='Turkana')
pval2<-c()
for (i in 1:10){
pval2<-c(pval2,wilcox.test(tmp[,3+i]~tmp$V5.y)$p.value)}

plot(1:10,p.adjust(pval1,method='BH'),col='steelblue',pch=20,bty='n',xlab='PC',ylab='FDR corrected p-value',cex.lab=1.2,cex.axis=1.2,ylim=c(0,1))
points(1:10,p.adjust(pval2,method='BH'),pch=20,col='goldenrod')
abline(h=0.1,lty=2)

