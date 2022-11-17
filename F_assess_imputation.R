##########
# RUN IN BASH
# concat imputation data for one chromosome
##########

touch /Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses/impute1panel_chr1.hg19_GQ10_info

for f in /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.*.hg19_GQ10_info; do awk '{OFS="\t"; print $0,FILENAME}' $f | tail -n+2 >> /Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses/impute1panel_chr1.hg19_GQ10_info; done

touch /Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses/impute2panel_chr1.hg19_GQ10_info

for f in /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute2panel.1.*.hg19_GQ10_info; do awk '{OFS="\t"; print $0,FILENAME}' $f | tail -n+2 >> /Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses/impute2panel_chr1.hg19_GQ10_info; done

cat /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.*.hg19_GQ10_info_by_sample | grep -v 'concord_type'> /Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses/impute1panel_chr1.hg19_GQ10_info_by_sample

cat /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute2panel.1.*.hg19_GQ10_info_by_sample | grep -v 'concord_type'> /Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses/impute2panel_chr1.hg19_GQ10_info_by_sample

# sample info
# cp /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/chr1.hg19.low_cov.GQ10.gen.samples /Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses

##########
# RUN IN R
# check imputation data for one chromosome
##########

setwd("~/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/Figures")

other_groups=c('1022','1035','1036','1038','1041','1044_Pokot','1046','1047','1048','1051','1052','1055','1056','1058','1060','1063','1070','1090','1091','1095','1098','1101','1103','1105','1106','1170','1182','1185','1196','1198','1212','1220','1224','1242','1244','1263','1266','1267','1269','1274','1275','1282','1296','1297','1299','1310','1311','1312','596','598','599','603','608','609','613','614','971','974','980','981_Samburu','985','987','997')

het=read.delim('~/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data/chr1.low_cov.GQ10.het.txt')
miss=read.delim('~/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data/chr1.low_cov.GQ10.imiss.txt')

keep<-subset(het,F> -0.2 & F< 0.2)
keep$id<-paste(keep$FID, keep$IID,sep='_')
miss$id<-paste(miss$FID, miss$IID,sep='_')
keep=merge(keep,miss,by='id')

# snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0
data1=read.delim('impute1panel_chr1.hg19_GQ10_info',sep=' ',header=F)
data2=read.delim('impute1panel_chr1.hg19_GQ10_info_by_sample',sep=' ',header=F)
# snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type1 concord_type1 r2_type1 info_type0 concord_type0 r2_type0
data3=read.delim('impute2panel_chr1.hg19_GQ10_info',sep=' ',header=F)
data4=read.delim('impute2panel_chr1.hg19_GQ10_info_by_sample',sep=' ',header=F)
data5=read.delim('chr1.hg19.low_cov.GQ10.gen.samples',sep=' ',header=F)
data5$order<-1:dim(data5)[1]

data2$order<-1:dim(data5)[1]
data4$order<-1:dim(data5)[1]

data5$conc_2panel=aggregate(data4[,3]~data4$order,FUN=median)[,2]
data5$r2_2panel=aggregate(data4[,4]~data4$order,FUN=median)[,2]

data5$conc_1panel=aggregate(data2[,1]~data2$order,FUN=median)[,2]
data5$r2_1panel=aggregate(data2[,2]~data2$order,FUN=median)[,2]

keep2=merge(data5,keep,by.y='id',by.x='V1',all.y=T)
keep2$group<-'steelblue'
keep2$group[which(keep2$FID.y %in% other_groups)]<-'goldenrod3'

keep2$group2<-'Turkana'
keep2$group2[which(keep2$FID.y %in% other_groups)]<-'Non Turkana'

# assess sample level quality for 1 vs 2 panel approach
par(mfrow=c(2,3))
plot(keep2$F_MISS,keep2$conc_2panel,col=(keep2$group),ylim=c(0.92,1),xlab='Proportion of data to impute',ylab='Concordant genotypes, 2 panel approach',pch=20,bty='n',cex.lab=1.2,cex.axis=1.2)
plot(keep2$F_MISS,keep2$conc_1panel,col=(keep2$group),ylim=c(0.92,1),xlab='Proportion of data to impute',ylab='Concordant genotypes, 1 panel approach',pch=20,bty='n',cex.lab=1.2,cex.axis=1.2)
plot(keep2$conc_1panel,keep2$conc_2panel,col=(keep2$group),xlab='Concordant genotypes, 1 panel approach',ylab='Concordant genotypes, 2 panel approach',ylim=c(0.92,1),xlim=c(0.92,1),pch=20,bty='n',cex.lab=1.2,cex.axis=1.2)
x=c(0,1); y=c(0,1); abline(lm(y~x),lty=2)

wilcox.test(keep2$conc_1panel,keep2$conc_2panel)
wilcox.test(keep2$r2_1panel,keep2$r2_2panel)

library(ggplot2)
# ggplot( data=keep2,aes(x=group2, y=conc_2panel)) + geom_violin(width=1.4) + geom_boxplot(width=0.1, color="grey", alpha=0.2) +theme_bw(13)

# check for differences between turkana and non turkana, controlling for coverage
summary(lm(keep2$conc_1panel~keep2$F_MISS + keep2$group))
summary(lm(keep2$r2_1panel~keep2$F_MISS + keep2$group))


# assess SNP level quality

both=merge(data1,data3,by='V3',bty='n')
plot(density(data1[,11]),xlab='SNP R2, 1 panel approach',xlim=c(0,1),cex.lab=1.2,cex.axis=1.2,bty='n',main='')
plot(density(data1[,7]),xlab='SNP info metric, 1 panel approach',xlim=c(0,1),cex.lab=1.2,cex.axis=1.2,bty='n',main='')

