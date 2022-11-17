#!/usr/bin/env Rscript

# NOTE: CHROMNUM was replaced with the chromosome number

# change Plink frequency files to format needed for BayEnv

library(data.table)
data=fread('/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.CHROMNUM.merged_LDfilt2_w1000G.frq.strat')

output1<-matrix(ncol=8,nrow=dim(data)[1]/8)
output2<-matrix(ncol=8,nrow=dim(data)[1]/8)
loci<-unique(data$SNP)

library(reshape2)
data_wide <- dcast(data[,c('SNP','CLST','MAC')], SNP ~ CLST, value.var="MAC")
data_wide2 <- dcast(data[,c('SNP','CLST','NCHROBS')], SNP ~ CLST, value.var="NCHROBS")
data_wide2b=subset(data_wide2,YRI>0 & LWK>0 &Turkana>0 )
data_wide=subset(data_wide,SNP %in% data_wide2b$SNP)
identical(data_wide2b$SNP,data_wide$SNP)

data_wide3<-data_wide2b
data_wide3$El_Molo<-data_wide2b$El_Molo-data_wide$El_Molo
data_wide3$Ik<-data_wide2b$Ik-data_wide$Ik
data_wide3$LWK<-data_wide2b$LWK-data_wide$LWK
data_wide3[,5]<-data_wide2b[,5]-data_wide[,5]
data_wide3$Rendille<-data_wide2b$Rendille-data_wide$Rendille
data_wide3$Samburu<-data_wide2b$Samburu-data_wide$Samburu
data_wide3$Turkana<-data_wide2b$Turkana-data_wide$Turkana
data_wide3$YRI<-data_wide2b$YRI-data_wide$YRI

output3<-matrix(ncol=8,nrow=dim(data_wide3)[1]*2)

evenseq<-seq(1, dim(data_wide3)[1]*2, 2)
oddseq<-seq(2, dim(data_wide3)[1]*2, 2)

for (f in 1:length(evenseq)){
output3[evenseq[f],]<-as.matrix(data_wide3[f,-c(1,5)])
output3[oddseq[f],]<-as.matrix(data_wide[f,-c(1,5)])
}

write.table(output3,'28May20_BayEnv_chrCHROMNUM_wLWKYRI.txt',row.names=F,col.names=F,sep='\t',quote=F)
write.table(data_wide3$SNP,'28May20_BayEnv_SNPINFO_chrCHROMNUM_wLWKYRI.txt',row.names=F,col.names=F,sep='\t',quote=F)

