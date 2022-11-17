#!/usr/bin/env Rscript

library(rehh)
library(data.table)
library(vcfR)

for (i in 1:22) {
	hh <- data2haplohh(hap_file = paste("high_cov.hg19_chr",i,".phased_withAA.vcf.gz",sep=''),polarize_vcf = TRUE,min_perc_geno.mrk=90)
	res<-scan_hh(hh,threads=4,phased = TRUE)
	res.ihs<-ihh2ihs(res)
	res.ihs.df2<-as.data.frame(res.ihs$ihs)
	res.ihs.df2<-res.ihs.df2[complete.cases(res.ihs.df2),]
	write.table(res.ihs.df2,paste("26May20_iHS_test_chr",i,".txt",sep=''),row.names=F,col.names=F,sep='\t',quote=F)
 	print(i) }