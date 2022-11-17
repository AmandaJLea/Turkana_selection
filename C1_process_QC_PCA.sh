###########
# RUN IN BASH
###########

# filter merged data for all array batches
~/programs/plink_1.90 --bfile ~/Turkana/genotyping_Jul22/merge_9batches --make-bed --out ~/Turkana/genotyping_Jul22/merge_9batches_filt1 --maf 0.01 --mind 0.05 --geno 0.05 --hwe 0.000001 --chr 1-22 --indep-pairwise 50 20 0.8 

# check relatedness
~/programs/king -b ~/Turkana/genotyping_Jul22/merge_9batches_filt1.bed --related --degree 2
~/programs/king -b ~/Turkana/genotyping_Jul22/merge_9batches_filt1.bed --unrelated --degree 2

# check het
~/programs/plink_1.90 --bfile ~/Turkana/genotyping_Jul22/merge_9batches_filt1 --het  --extract ~/Turkana/genotyping_Jul22/merge_9batches_filt1.prune.in

############
# RUN IN R
############

# check dups
data=read.delim('king.kin0',stringsAsFactors=FALSE)
data2=subset(data,InfType=='Dup/MZ')
write.table(data2,'duplicates.txt',row.names=F,sep='\t')

# prep for removal
remove<-c('273','2073','2089_4','370','2169','1604_5','1853_7','1882','1766_9','1894','786','776_7','2240','272','1642_6','1627','1354_7','1356_9','2269','308','1563','1865','1896','844_11','931_1','969_10')

data=read.delim('kingunrelated_toberemoved2.txt',header=F)
fam=read.delim('~/Turkana/genotyping_Jul22/merge_9batches_filt1.fam',header=F,sep=' ')
fam2=subset(fam,(V2 %in% data$V1) | V2 %in% remove)
write.table(fam2,'kingunrelated_toberemoved3.txt',row.names=F,sep='\t',col.names=F,quote=F)

############
# RUN IN BASH
############

~/programs/plink_1.90 --bfile ~/Turkana/genotyping_Jul22/merge_9batches_filt1 --out ~/Turkana/genotyping_Jul22/merge_9batches_filt1 --remove kingunrelated_toberemoved3.txt --maf 0.01 --pca 10 tabs --extract ~/Turkana/genotyping_Jul22/merge_9batches_filt1.prune.in

############
# RUN IN R
############

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures")

meta=read.delim('15Jul22_all_genotyped_samples.txt')
data1=read.delim('merge_9batches_filt1.eigenval',header=F)
data2=read.delim('merge_9batches_filt1.eigenvec',header=F)

data3=merge(meta,data2,by.x='FileID',by.y='V2')
write.table(data3,'18Aug22_array_data_filtered.txt',row.names=F,sep='\t')

library(ggplot2)
ggplot(data3, aes(x=V3, y=V4, col=Group)) +geom_point(alpha=0.65,size=2)+theme_bw(15)+xlab("PC 1")+ylab("PC 2")+ scale_color_brewer(palette="Paired")

plot(1:10,data1$V1,xlab='PC',ylab='PVE',bty='n',pch=20)

ggplot(data3, aes(x=Group, y=V4, col=Group)) +geom_violin(width=2)+geom_boxplot(width=0.4,, alpha=0.2)+theme_bw(15)+xlab("Group")+ylab("PC 1")+ scale_color_brewer(palette="Paired")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# R2
x1<-apply(data3[,c(7:16)],2,function(x) summary(lm(x~data3$Group))$r.squared)

# sampling locations
pheno=read.delim("18Aug22_locations.txt")
data4=merge(pheno,data3,by.x='Individual.sbarcodenumber.x',by.y='FileID',all.y=T)

# el molo
data4$X.y[which(data4$Group=='El Molo')]<-36.41530
data4$Y[which(data4$Group=='El Molo')]<-2.51328

data4=subset(data4, X.y> -100 & Group=='Turkana')

tmp<-data4[,c('X.y','Y','V3','V4','V5','V6','V7')]
names(tmp)<-c('Longitude','Latitude','PC1','PC2','PC3','PC4','PC5')
library(corrplot)
corrplot(cor(tmp),method='ellipse',type='upper')

ggplot(tmp, aes(x=Longitude, y=PC1)) +geom_point(alpha=0.65,size=2)+theme_bw(15)+xlab("Longitude")+ylab("PC 1")+ geom_smooth(method = "lm") 
ggplot(tmp, aes(x=Latitude, y=PC2)) +geom_point(alpha=0.65,size=2)+theme_bw(15)+xlab("Latitude")+ylab("PC 2")+ geom_smooth(method = "lm") 

# R2
tmp<-data4[,c('X.y','Y','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')]

x2<-apply(tmp[,c(3:12)],2,function(x) summary(lm(x~tmp$X.y + tmp$Y))$r.squared)

par(mfrow=c(2,2))
plot(1:10,x1,pch=20,col='steelblue',xlab='PC',bty='n',ylab='PVE by group',ylim=c(0,0.45))
plot(1:10,x2,pch=20,col='steelblue',xlab='PC',bty='n',ylab='PVE by location',ylim=c(0,0.45))
