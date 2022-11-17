#!/usr/bin/env Rscript

# NOTE: CHROMNUM was replaced with the chromosome number

# calculate PBS

library(data.table)

x<-fread("YRI_Turk_chrCHROMNUM_FST.fst",header=T)
y<-fread("LWK_Turk_chrCHROMNUM_FST.fst",header=T)
z<-fread("LWK_YRI_chrCHROMNUM_FST.fst",header=T)

identical(x$SNP,y$SNP)
identical(z$SNP,y$SNP)

x$T1<- -log(1-x$FST) 
x$T2<- -log(1-y$FST) 
x$T3<- -log(1-z$FST) 

x$pbs_lwk_yri<- (x$T2 + x$T1 - x$T3)/2

tmp<-read.delim('LWK_YRI_Turk_chrCHROMNUM.frq_v2.strat')
library(reshape2)
data_wide <- dcast(tmp, SNP ~ CLST, value.var="MAF")

out<-merge(x,data_wide,by='SNP')
write.table(out[,-c(5:8)],'19May20_PBS_test_chrCHROMNUM.txt',row.names=F,sep='\t',quote=F)

