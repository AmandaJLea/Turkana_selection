
##########
# combined selection scores into a Fisher’s score (FCS) equal to the sum, over the two statistics, of –log10(rank of the statistic for a given SNP/number of SNPs)
##########

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

all=unique(merge(res1,res2[,c('snp','iHS_stat')],by='snp'))
all=unique(merge(all,data1[,c('V1','mean_stat')],by.x='snp',by.y='V1'))
all$abs_iHS_stat<-abs(all$iHS_stat)
all=all[complete.cases(all),]
all=unique(all[,-5])

all=(all[order(all$abs_iHS_stat,decreasing =T),])
all$abs_iHS_stat_rank<-1:dim(all)[1]

all=(all[order(all$pbs_stat,decreasing =T),])
all$pbs_stat_rank<-1:dim(all)[1]

all=(all[order(all$mean_stat,decreasing =T),])
all$xtx_stat_rank<-1:dim(all)[1]

all$fcs<- -log10(all$abs_iHS_stat_rank/1578201)+ -log10(all$pbs_stat_rank/1578201) + -log10(all$xtx_stat_rank/1578201)
all$fcs2<- -log10(all$abs_iHS_stat_rank/1578201) + -log10(all$xtx_stat_rank/1578201)
all$loc2<-all$V2+1

write.table(all,'/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/2Aug20_combined_selection_FCS.txt',row.names=F,sep='\t')
write.table(all[,c('V1','V2','loc2','fcs2')],'/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/2Aug20_combined_selection_FCS.bed',row.names=F,sep='\t',quote=F,col.names=F)

##########
# computed for each genomic window, associated or not with the trait, the average FCS, the proportion of conserved SNP positions based on GERP scores > 2 [92], and the recombination rate using the combined HapMap genetic map [102], to account for the confounding effects of background selection.
##########

rm /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements.bed; touch /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements.bed

for f in {1..22}; do awk -v CHR=$f '{OFS=""; print CHR,"\t",$1,"\t",$2}' /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_chr${f}_elems.txt >> /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements.bed; done

bedtools intersect -a /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_all_elements.bed -b $windows -wo > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_50kb_windows_overlap.txt

rm /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_all_chr.txt; touch /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_all_chr.txt

for f in {1..22}; do awk -v CHR=$f '{OFS=""; print CHR,"\t",$1,"\t",$1+1,"\t",$2}' /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_chr${f}_combined_b37.txt | tail -n +2 >> /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_all_chr.txt; done

bedtools intersect -a /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/genetic_map_all_chr.txt -b $windows -wo > /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/hg19_50kb_windows_overlap.txt

bedtools intersect -a /Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/2Aug20_combined_selection_FCS.bed -b $windows -wo > /Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/2Aug20_combined_selection_FCS_50kb_windows_overlap.bed
#

# non overlapping windows to keep
no_overlap<-read.delim('/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_v0_Homo_sapiens_assembly19.50kb_windows.nooverlap.txt',header=F)
no_overlap$V2[which(no_overlap$V2==0)]<-1
no_overlap$window<-paste(no_overlap$V1,no_overlap$V2,no_overlap$V3,sep='_')

library(data.table)
recomb=fread('/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000GP_Phase3/hg19_50kb_windows_overlap.txt',header=F)
gerp=fread('/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/GERP_conserved_elements/hg19_50kb_windows_overlap.txt',header=F)
fcs=fread('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/2Aug20_combined_selection_FCS_50kb_windows_overlap.bed',header=F)

fcs$window<-paste(fcs$V5,fcs$V6,fcs$V7,sep='_')
gerp$window<-paste(gerp$V4,gerp$V5,gerp$V6,sep='_')
recomb$window<-paste(recomb$V5,recomb$V6,recomb$V7,sep='_')

# fcs by window
tmp1<-aggregate(fcs$V4~fcs$window,FUN=mean)
tmp2<-aggregate(fcs$V4~fcs$window,FUN=length)
identical(tmp1[,1],tmp2[,1])
fcs_sum<-as.data.frame(cbind(tmp1[,2],tmp2[,2]))
fcs_sum$window<-tmp1[,1]

# conserved bases by window
gerp_sum<-aggregate(gerp$V7~gerp$window,FUN=sum)
names(gerp_sum)<-c('window','conserved_bp')

# recomb by window
tmp1<-aggregate(recomb$V4~recomb$window,FUN=sum)
tmp2<-aggregate(recomb$V4~recomb$window,FUN=length)
identical(tmp1[,1],tmp2[,1])
recomb_sum<-as.data.frame(cbind(tmp1[,2],tmp2[,2]))
recomb_sum$window<-tmp1[,1]
recomb_sum$mean_recomb<- recomb_sum$V1 /50000

all=merge(recomb_sum,gerp_sum,by='window',all.x=T,all.y=T)
all=merge(all,fcs_sum,by='window',all.x=T,all.y=T)
names(all)<-c("window","recomb_sum","recomb_count","mean_recomb","conserved_bp","mean_fcs","fcs_snps") 

all2=all[which(all$window %in% no_overlap$window),]

write.table(all2,'/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/2Aug20_combined_selection_FCS_wGERP_wrecomb_nooverlap.txt',row.names=F,sep='\t')

##########
# to test for polygenic selection, we generated a null distribution by randomly sampling x windows (x being the number of windows associated with a tested trait) among windows with a similar number of SNPs, proportion of GERP > 2 sites and recombination rate observed in the trait-associated windows. 
# We then calculated the average of the mean of the FCS across the x resampled windows. 
# We resampled 100,000 sets of x windows for each trait. 
##########

# cd /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/UKBB_pan
# for f in windows.trait*sig.bed ; do cat temp.sh | sed -e s/FILEINFO/${f}/g > temp.${f}.R; rm temp.${f}.sh; touch temp.${f}.sh; echo '#!/bin/bash' >> temp.${f}.sh; echo "module load R/3.5.2; Rscript temp.${f}.R" >> temp.${f}.sh; done
# ls *-sig.bed.sh | sed -e 's/temp.windows.trait/sh temp.windows.trait/g' > commands1.sh

fcs=read.delim('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/2Aug20_combined_selection_FCS_wGERP_wrecomb_nooverlap.txt',header=T)
fcs=subset(fcs,fcs_snps>0)

tmp<-quantile(fcs$mean_recomb,na.rm=T)
fcs$recomb_quant<-NA
fcs$recomb_quant[which(fcs$mean_recomb <= tmp[5])]<-4
fcs$recomb_quant[which(fcs$mean_recomb <= tmp[4])]<-3
fcs$recomb_quant[which(fcs$mean_recomb <= tmp[3])]<-2
fcs$recomb_quant[which(fcs$mean_recomb <= tmp[2])]<-1

tmp<-quantile(fcs$conserved_bp,na.rm=T)
fcs$gerp_quant<-NA
fcs$gerp_quant[which(fcs$conserved_bp <= tmp[5])]<-4
fcs$gerp_quant[which(fcs$conserved_bp <= tmp[4])]<-3
fcs$gerp_quant[which(fcs$conserved_bp <= tmp[3])]<-2
fcs$gerp_quant[which(fcs$conserved_bp <= tmp[2])]<-1

tmp<-quantile(fcs$fcs_snps,na.rm=T)
fcs$num_quant<-NA
fcs$num_quant[which(fcs$fcs_snps <= tmp[5])]<-4
fcs$num_quant[which(fcs$fcs_snps <= tmp[4])]<-3
fcs$num_quant[which(fcs$fcs_snps <= tmp[3])]<-2
fcs$num_quant[which(fcs$fcs_snps <= tmp[2])]<-1

sig=read.delim('FILEINFO',header=F)
sig$window<-paste(sig$V1,sig$V2,sig$V3,sep='_')

sig_info<-fcs[which(fcs$window %in% sig$window),]
sig_info2<-as.data.frame(table(sig_info$recomb_quant,sig_info$gerp_quant,sig_info$num_quant))
notsig_info<-fcs[-which(fcs$window %in% sig$window),]

null_sample<-c()

# null distribution
output<-lapply(c(1:1000),function(x) {

null_sample<-lapply(1:dim(sig_info2)[1], function(i) {
tmp1<-(subset(notsig_info,recomb_quant==sig_info2$Var1[i] & gerp_quant==sig_info2$Var1[i] & num_quant==sig_info2$Var1[i]))
return(tmp1[sample(1:dim(tmp1)[1],sig_info2$Freq[i]),]) })

mean_tmp<-c()
for (k in 1:dim(sig_info2)[1]){
mean_tmp<-c(mean_tmp,null_sample[[k]][,6])
}
return(mean(mean_tmp))
} )

write.table(as.matrix(output),'results2_FILEINFO',row.names=F,sep='\t')

##########
# To test for significance, we computed a resampling P-value by calculating the proportion of resampled windows which mean FCS was higher than that observed for the tested trait. All P-values for polygenic adaptation were then adjusted for multiple testing by the Benjamini-Hochberg method, to account for the number of traits tested, and traits with an adjusted p < 0.05 were considered as candidates for polygenic selection.
##########

# cd /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/UKBB_pan

files2=read.delim('files.txt',header=F)

fcs=read.delim('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/2Aug20_combined_selection_FCS_wGERP_wrecomb_nooverlap.txt',header=T)
fcs=subset(fcs,fcs_snps>0)

remove<-c('12_130175000_130225000','12_130200000_130250000','18_66700000_66750000','3_21150000_21200000','5_115850000_115900000','5_115875000_115925000','5_115900000_115950000','8_15075000_15125000','8_15100000_15150000','8_23850000_23900000','8_23875000_23925000','8_3125000_3175000','9_1500000_1550000')

fcs=fcs[-which(fcs$window %in% remove),]

pval<-c()
mean_null<-c()
mean_obs<-c()
sd_null<-c()

for (i in c(1:22,24:30)){
x=paste("results2_windows.trait-",files2[i,1],"-sig.bed",sep="")
y=read.delim(x)

sig=read.delim(paste('windows.trait-',files2[i,1],'-sig.bed',sep=''),header=F)
sig$window<-paste(sig$V1,sig$V2,sig$V3,sep='_')

sig2<-fcs[which(fcs$window %in% sig$window),]

pval<-c(pval,length(which(y$V1 < mean(sig2$mean_fcs))) /1000)
mean_null<-c(mean_null,mean(y$V1))
mean_obs<-c(mean_obs,mean(sig2$mean_fcs))
sd_null<-c(sd_null,sd(y$V1))

}

output<-cbind(mean_null,sd_null,mean_obs,1-pval,p.adjust(1-pval,method='BH'))



