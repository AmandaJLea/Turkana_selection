############
# furcation plot
############

setwd('/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures')

i<-c('high_cov_chr8_STC1.vcf')

# 8:23872028, 8:23873255, 8:23873475, 8:23875208

library(rehh)

	hh <- data2haplohh(hap_file = i,polarize_vcf = TRUE,min_perc_geno.mrk=90)
	res<-scan_hh(hh,threads=4,phased = TRUE)
	res.ihs<-ihh2ihs(res)
	res.ihs.df2<-as.data.frame(res.ihs$ihs)
	res.ihs.df2<-res.ihs.df2[complete.cases(res.ihs.df2),]

	
	lead<-res.ihs.df2[which(res.ihs.df2$POSITION==23872028),]
	furc <- calc_furcation(hh, mrk = which(hh@positions == lead$POSITION))
	# plot(furc, xlim = c(lead$POSITION-20000 , lead$POSITION+10000),main='',legend = NA )
	
haplen <- calc_haplen(furc)
plot(haplen,main=" ",legend = NA,col=c('darkseagreen','steelblue') )

########
# mouse ATAC-seq
########

# hg19 - region
# 8	23850000	23925000

# hg 19 - STC1
# 8 23699428	23712320

# mm10 - region
# chr14	68860352	68921736

# mm10 - STC1
# chr14	69029125	69041405

setwd('/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures')

data1=read.delim('GSE108786_RAW/GSM2913206_Rep1_vehicle_ATACSeq_macs2_peak.bed',header=F)
data2=read.delim('GSE108786_RAW/GSM2913207_Rep1_dDAVP_ATACSeq_macs2_peak.bed',header=F)
data3=read.delim('GSE108786_RAW/GSM2913208_Rep2_vehicle_ATACSeq_macs2_peak.bed',header=F)
data4=read.delim('GSE108786_RAW/GSM2913209_Rep2_dDAVP_ATACSeq_macs2_peak.bed',header=F)

data1$id<-'-dDAVP'
data2$id<-'+dDAVP'
data3$id<-'-dDAVP'
data4$id<-'+dDAVP'

all<-rbind(data1,data2,data3,data4)

count<-c()
selected<-read.delim('selected_regions_mm10.txt',header=F)
selected$length<-selected$V3-selected$V2
for (i in 1:dim(selected)[1]){
tmp<-subset(all,V1==selected$V1[i] & V2>selected$V2[i] & V3<selected$V3[i] )
count<-c(count,dim(tmp)[1])
}

all2A<-subset(all,V1=='chr14' & V2>68860352-10000 & V3<69041405+10000 & id=='-dDAVP')
all2B<-subset(all,V1=='chr14' & V2>68860352-10000 & V3<69041405+10000 & id=='+dDAVP')

coords<-c(68860352-10000,68921736,69029125,69041405+10000)
plot(1, type = "n", xlab = "", ylab = "",xlim=c(min(coords),max(coords)),ylim=c(0.8,3),bty='n')

 segments(min(coords), 0.9, max(coords), 0.9,lwd=2,col='grey90')
 segments(min(coords), 2, max(coords), 2,lwd=2,col='grey90')

segments(69029125, 0.9, 69041405, 0.9,lwd=5,col='darkseagreen')
segments(69029125, 0.9, 69029125, 1.3,lwd=5,col='darkseagreen')
segments(69029125, 1.3, 69029125+1000, 1.3,lwd=5,col='darkseagreen')
arrows(69029125+8000, 1.3, 69029125+8100, 1.3,lwd=5,col='darkseagreen',length = 0.1)
segments(69029125, 1.3, 69029125+8100, 1.3,lwd=5,col='darkseagreen')

 segments(68860352, 0.9, 68921736, 0.9,lwd=2,col='darkgrey')

segments(69029125, 2, 69041405, 2,lwd=5,col='darkseagreen')
segments(69029125, 2, 69029125, 2.4,lwd=5,col='darkseagreen')
segments(69029125, 2.4, 69029125+1000, 2.4,lwd=5,col='darkseagreen')
arrows(69029125+8000, 2.4, 69029125+8100, 2.4,lwd=5,col='darkseagreen',length = 0.1)
segments(69029125, 2.4, 69029125+8100, 2.4,lwd=5,col='darkseagreen')

segments(68860352, 2, 68921736, 2,lwd=2,col='darkgrey')
# abline(v=68860352, lty=2,col='grey')
# abline(v=68921736, lty=2,col='grey')

for (i in 1:dim(all2A)[1]){
segments(all2A$V2[i], 1.1, all2A$V3[i], 1.1,lwd=5,col='skyblue')
}

for (i in 1:dim(all2B)[1]){
segments(all2B$V2[i], 2.2, all2B$V3[i], 2.2,lwd=5,col='steelblue3')
}

########
# qPCR
########

setwd('/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures')

data=read.delim('20Sep22_qPCR_data.txt')
data=subset(data,Use=='Y')
data$Label3 <- factor(data$Label2, levels=c('PBS','10nM','20nM'))

treat<-subset(data,Label2!='PBS')
ctrl<-subset(data,Label2=='PBS')
x<-mean(ctrl$Delta_Ct)
treat$Fold_change<-2^-(treat$Delta_Ct - x)

library(ggplot2)
ggplot(treat, aes(x=Label2, y=Fold_change)) + geom_boxplot()+geom_jitter()+theme_bw(15)+ylab('Fold change in expression')+coord_flip()

ggplot(data, aes(x=Label3, y=Delta_Ct)) + geom_boxplot()+geom_jitter()+theme_bw(15)+ylab('Delta Ct')+coord_flip()


t.test(subset(treat,Label3=='10nM')$Delta_Ct,ctrl$Delta_Ct)
t.test(subset(treat,Label3=='20nM')$Delta_Ct,ctrl$Delta_Ct)