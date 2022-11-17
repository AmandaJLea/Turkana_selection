##########
# combined selection scores
##########

setwd('/Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses')
i=1
data1=read.delim(paste('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/bayenv_files/8Jun20_BayEnv_results_chr',i,'.txt',sep=''))
data2=read.delim(paste('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/bayenv_files/8Jun20_BayEnv_region_results_chr',i,'.txt',sep=''))

for (i in 2:22){
tmp1=read.delim(paste('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/bayenv_files/8Jun20_BayEnv_results_chr',i,'.txt',sep=''))
tmp2=read.delim(paste('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/bayenv_files/8Jun20_BayEnv_region_results_chr',i,'.txt',sep=''))
names(tmp1)<-names(data1)
names(tmp2)<-names(data2)
data1=rbind(data1,tmp1)
data2=rbind(data2,tmp2)
print(i)
}

library(data.table)
res1=fread('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/PBS_files/19May20_PBS_test_ALLchr.txt',header=F)
res2=fread('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/rehh_files/26May20_iHS_test_ALLchr.bed',header=F)
res1$snp<-paste(res1$V1,res1$V2,sep=':')
res2$snp<-paste(res2$V1,res2$V2,sep=':')
names(res1)[3]<-'pbs_stat'
names(res2)[3]<-'iHS_stat'

# make bed file of all 3 stats

by_snp<-merge(res1,data1[,c('V1','mean_stat')],by.x='snp',by.y='V1',all.x=T,all.y=T)
by_snp2<-merge(by_snp,res2[,-c(1:2)],by='snp',all.x=T,all.y=T)
write.table(by_snp2[,c('snp','iHS_stat','mean_stat','pbs_stat')],'24Feb22_3selection_stats_by_SNP.txt',row.names=F,sep='\t',quote=F)

##

cat 24Feb22_3selection_stats_by_SNP.txt | sed -e s/:/"\t"/g | tail -n+2 | awk '{OFS="\t"; print $1,$2,$2,$3,$4,$5}' > 24Feb22_3selection_stats_by_SNP.bed

# note some duplicated rows

##########
# link up with genes - hg19
##########

# add 100kb to protein coding
cat /Genomics/ayroleslab2/alea/ref_genomes/hg19/Homo_sapiens.GRCh37.87_protein_coding_sort.bed | awk '{OFS="\t"; print $1,$2-100000,$3+100000,$5}' | awk '$2>1' | awk '$3>1' | sed -e 's/"//g' | sed -e 's/;//g' > /Genomics/ayroleslab2/alea/ref_genomes/hg19/Homo_sapiens.GRCh37.87_protein_coding_100kb.bed

module load bedtools
intersectBed -a /Genomics/ayroleslab2/alea/ref_genomes/hg19/Homo_sapiens.GRCh37.87_protein_coding_100kb.bed -b 24Feb22_3selection_stats_by_SNP.bed -wo >24Feb22_3selection_stats_by_SNP_wgene.bed

##########
# link up with DE genes
##########

cand<-c('ENSG00000151789','ENSG00000092421','ENSG00000183117','ENSG00000185053','ENSG00000159167','ENSG00000173253','ENSG00000151952','ENSG00000150636')

setwd("/Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses")
library(data.table)
selection=fread('24Feb22_3selection_stats_by_SNP_wgene.bed',header=F)
selection=unique(selection)
exp=read.delim('20Sep22_LMM_pvals.txt')
exp2=subset(exp,env_qvalue<0.1)
exp1=subset(exp,env_qvalue>0.2)

# exp level per gene
data2=fread('/Genomics/ayroleslab2/alea/turkana_wgs/24Aug21_feature_counts.txt')
tot_counts<-apply(data2[,-1],2,sum)
remove<-which(tot_counts<250000)
data3<-data2[,-c(remove+1),with=F]
tot_counts<-apply(data3[,-1],2,sum)
norm1=10^6 *(data3[,-1]/tot_counts)
data3$med_TPM<-apply(norm1,1,median)
data3$mean_TPM<-apply(norm1,1,mean)
data3=subset(data3,Gene %in% exp$geneID)

n=10
qtile = quantile(data3$med_TPM, probs = seq(0, 1, 1/n))
data3$exp_quartile = sapply(data3$med_TPM, function(x) sum(x >= qtile[-(n+1)]))

selection2=merge(selection,data3[,c('Gene','med_TPM','exp_quartile')],by.x='V4',by.y='Gene')
selection2$DE<-NA
selection2$DE[which(selection2$V4 %in% exp2$gene)]<-1
selection2$DE[which(selection2$V4 %in% exp1$gene)]<-0

# per SNP - near selected vs un selected
tmp<-unique(selection2[,c('V8','DE','exp_quartile','V6')])
summary(lm( abs(tmp$V8)~tmp$DE+tmp$exp_quartile))
x1<-ggplot(tmp, aes(x=V8, color=factor(DE))) +geom_density()+scale_fill_brewer(palette='Paired')

tmp<-unique(selection2[,c('V9','DE','exp_quartile','V6')])
summary(lm(tmp$V9~tmp$DE+tmp$exp_quartile))
X2<-ggplot(tmp, aes(x=V8, color=factor(DE))) +geom_density()+scale_fill_brewer(palette='Paired')

tmp<-unique(selection2[,c('V10','DE','exp_quartile','V6')])
summary(lm(tmp$V10~tmp$DE+tmp$exp_quartile))
ggplot(tmp, aes(x=V8, color=factor(DE))) +geom_density()+scale_fill_brewer(palette='Paired')
       
# per gene
tmp<-unique(selection2[,c('V8','DE','exp_quartile','V6','V4')])
selection3<-as.data.frame(aggregate( abs(tmp$V8)~tmp$V4,FUN=median))
tmp<-unique(selection2[,c('V9','V10','DE','exp_quartile','V6','V4')])
selection4<-as.data.frame(aggregate(tmp$V9~tmp$V4,FUN=median))
selection4$xtx<-as.data.frame(aggregate(tmp$V10~tmp$V4,FUN=median))[,2]
names(selection3)<-c('gene','ihs')
names(selection4)<-c('gene','pbs','xtx')
selection5<-merge(selection3,selection4,by='gene')
selection5<-merge(selection5,data3[,c('Gene','exp_quartile')],by.x='gene',by.y='Gene')
selection5$DE<-NA
selection5$DE[which(selection5$gene %in% exp2$gene)]<-1
selection5$DE[which(selection5$gene %in% exp1$gene)]<-0

summary(lm(selection5$ihs~selection5$DE+selection5$exp_quartile))
summary(lm(selection5$pbs~selection5$DE+selection5$exp_quartile))
summary(lm(selection5$xtx~selection5$DE+selection5$exp_quartile))

# outliers
selection2$ihs_outlier<-0
selection2$ihs_outlier[which( abs(selection2$V8) > quantile(abs(selection2$V8),seq(0,1,0.01),na.rm=T)[100])]<-1

selection2$pbs_outlier<-0
selection2$pbs_outlier[which( (selection2$V9) > quantile((selection2$V9),seq(0,1,0.01),na.rm=T)[100])]<-1

selection2$xtx_outlier<-0
selection2$xtx_outlier[which( (selection2$V10) > quantile((selection2$V10),seq(0,1,0.01),na.rm=T)[100])]<-1

selection3<-as.data.frame(aggregate( abs(selection2$ihs_outlier)~selection2$V4,FUN=sum))
selection4<-as.data.frame(aggregate(selection2$pbs_outlier~selection2$V4,FUN=sum))
selection4$xtx<-as.data.frame(aggregate(selection2$xtx_outlier~selection2$V4,FUN=sum))[,2]
names(selection3)<-c('gene','ihs')
names(selection4)<-c('gene','pbs','xtx')
selection5<-merge(selection3,selection4,by='gene')
selection5<-merge(selection5,data3[,c('Gene','exp_quartile')],by.x='gene',by.y='Gene')
selection5$DE<-NA
selection5$DE[which(selection5$gene %in% exp2$gene)]<-1
selection5$DE[which(selection5$gene %in% exp1$gene)]<-0

summary(lm(selection5$ihs~selection5$DE+selection5$exp_quartile))
summary(lm(selection5$pbs~selection5$DE+selection5$exp_quartile))
summary(lm(selection5$xtx~selection5$DE+selection5$exp_quartile))
