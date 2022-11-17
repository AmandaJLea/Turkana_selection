##########
# RUN IN BASH
##########

# hg19 - region
# 8	23850000	23925000

# hg 19 - STC1
# 8 23699428	23712320

# temp.txt
# 8 23699428 23925000 1

cd /Genomics/ayroleslab2/alea/turkana_wgs/Jul2022_analyses

plink=/Genomics/grid/users/alea/programs/plink_1.90
path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets

# get AF in candidate region for 1KG
$plink --vcf $path_resources/shapeit_files/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --freq --out STC1_allele_freq_1KG --within /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/1000_genomes_pops.txt --extract range temp.txt

$plink --vcf $path_resources/shapeit_files/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-just-bim --out STC1_allele_freq_1KG --extract range temp.txt

# get AF in candidate region for Turkana
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.8.merged_LDfilt_plink --freq --out STC1_allele_freq_Turkana.txt --within turkana_data_groups.txt --extract range temp.txt

##########
# RUN IN R
##########

library(ggplot2)

# plot
setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures")
library(data.table)
data1=fread('STC1_allele_freq_1KG.frq.strat')
data2=fread('STC1_allele_freq_Turkana.txt.frq.strat.txt')
data1_bim=fread('STC1_allele_freq_1KG.bim')

data1<-merge(data1,data1_bim[,c(2,4),with=F],by.x='SNP',by.y='V2')
data1=subset(data1,A1 %in% c('A','T','G','A') & A2 %in% c('A','T','G','A'))
data2=subset(data2,A1 %in% c('A','T','G','A') & A2 %in% c('A','T','G','A'))

info=read.delim("1KG_population_codes.txt")
data1=merge(data1,info[,c(1,3)],by.x='CLST',by.y='Population.Code')
data1$SNP2<-paste(data1$CHR,data1$V4,sep=':')

both=merge(data1,data2,by.x='SNP2',by.y='SNP')

# MAF 
tmp<-subset(data1, SNP %in% both$SNP & V4>23850000	& V4<23925000 & (Super.Population.Code =='EUR' | Super.Population.Code=='AFR'))
dim(subset(tmp,MAF>0.1 & CLST=='CEU'))
ggplot(tmp, aes(x=V4, y=MAF,col=CLST)) + geom_point() + facet_wrap(~Super.Population.Code)+scale_color_brewer(palette='Paired')+theme_bw(13)
tmp<-subset(data2, SNP %in% both$SNP2  )
dim(subset(tmp,MAF>0.1 & CLST=='Turkana'))
ggplot(tmp, aes(x=LOC, y=MAF,col=CLST)) + geom_point() +scale_color_brewer(palette='Paired')+theme_bw(13)

both1=subset(both,A1.x==A1.y & A2.x==A2.y)
both2=subset(both,A1.x==A2.y & A1.y==A2.x)

data1_sub<-subset(data1,V4==23872028)
data2_sub<-subset(data2,SNP=='8:23872028')
data2_sub$Super.Population.Code<-'NWK'

both_sub<-rbind(data1_sub[,c('CLST','A1','A2','MAF','Super.Population.Code')],data2_sub[,c('CLST','A1','A2','MAF','Super.Population.Code')])
ggplot(both_sub, aes(x=Super.Population.Code, y=1-MAF)) + geom_boxplot()+geom_jitter( aes(x=Super.Population.Code, y=1-MAF,col=CLST))

# write.table(both_sub,'19Sep22_focal_SNP_MAF.txt',row.names=F,sep='\t')
data=read.delim('19Sep22_focal_SNP_MAF_wGPS.txt')

# plot maps
register_google('AIzaSyBqRWaHZTtsvaSOIwDdHjmTVPD-cXvaKrE')
library(ggplot2);library(ggmap);library(devtools);library(reshape2)

s <- "element:geometry%7Ccolor:0xf5f5f5&style=element:labels%7Cvisibility:off&style=element:labels.icon%7Cvisibility:off&style=element:labels.text.fill%7Ccolor:0x616161&style=element:labels.text.stroke%7Ccolor:0xf5f5f5&style=feature:administrative%7Celement:geometry%7Cvisibility:off&style=feature:administrative.country%7Celement:geometry.stroke%7Ccolor:0x000000%7Cvisibility:on&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.land_parcel%7Celement:labels.text.fill%7Ccolor:0xbdbdbd&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:poi%7Cvisibility:off&style=feature:poi%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:poi%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:poi.park%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:poi.park%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:geometry%7Ccolor:0xffffff&style=feature:road%7Celement:labels.icon%7Cvisibility:off&style=feature:road.arterial%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:road.highway%7Celement:geometry%7Ccolor:0xdadada&style=feature:road.highway%7Celement:labels.text.fill%7Ccolor:0x616161&style=feature:road.local%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:transit%7Cvisibility:off&style=feature:transit.line%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:transit.station%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:water%7Celement:geometry%7Ccolor:0xc9c9c9&style=feature:water%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&size=480x360"
map <- get_googlemap(center=c(-50,0), zoom = 2, scale = 1, style = s)
m<- ggmap(map)

library(scatterpie);library("ggrepel")

data$DAF<-1-data$MAF
data2=subset(data,(Long<20 & Super.Population.Code=='AFR') | CLST=='LWK')
m + geom_scatterpie(aes(x=Long, y=Lat), data=data2, alpha=0.75,cols=c( "DAF","MAF"),pie_scale=1.25) + scale_fill_manual(values=c('steelblue','goldenrod'))+ geom_label_repel(data=data2,aes(x=Long, y=Lat,label = CLST), size = 3.5)

map <- get_googlemap(center=c(36,2), zoom = 7, scale = 1, style = s)
m<- ggmap(map)
data2=subset(data,Long>20 & CLST!='LWK' )
data2[which(data2$CLST=='Masaai'),'Lat']<- -1
data2[which(data2$CLST=='Masaai'),'Long']<- 35
data2[which(data2$CLST=='Rendille'),'Lat']<- 2.2
data2[which(data2$CLST=='Rendille'),'Long']<- 36.75

m + geom_scatterpie(aes(x=Long, y=Lat), data=data2, alpha=0.75,cols=c( "DAF","MAF"),pie_scale=4) + scale_fill_manual(values=c('steelblue','goldenrod'))+ geom_label_repel(data=data2,aes(x=Long, y=Lat,label = CLST), size = 3.5,force=10)

data2=subset(data, CLST!='CEU' )
getPalette = colorRampPalette(brewer.pal(8, "Paired"))

ggplot(data2, aes(x=Super.Population.Code, y=1-MAF)) + geom_boxplot()+geom_jitter( aes(x=Super.Population.Code, y=1-MAF,col=CLST),size=2.5)+theme_bw(13)+scale_color_manual(values = getPalette(15))+ylab('DAF')

ggplot(data2, aes(x=Super.Population.Code, y=1-MAF)) + geom_boxplot()+geom_jitter( aes(x=Super.Population.Code, y=1-MAF),size=2.5)+theme_bw(15)+ylab('DAF')+coord_flip()


