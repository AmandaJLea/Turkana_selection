########
# ALL array data
########

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/array_genotyping/CHOP_batch9")

meta=read.delim('15Jul22_all_genotyped_samples.txt')
meta$type<-'array'
meta$team<-'Princeton'
meta_v2<-meta[,c('Barcode','Unique_ID..if.included.','Group','type','team')]

########
# ALL WGS data
########

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data")

data3=read.delim('impute1panel.1.merged_LDfilt_plink.fam_winfo.txt',header=F)
data3$V8<-'WGS'
data3_v2<-data3[,c('V1','V6','V3','V8','V4')]
names(data3_v2)<-names(meta_v2)
both<-rbind(meta_v2,data3_v2)

##########
# GPS locations + climate
##########

library(ggplot2);library(ggmap);library(devtools);library(reshape2)

# sampling locations
pheno=read.delim("~/Dropbox/Amanda_files/ayroles_lab/Turkana_project/early_life/data/9Oct20_cleaned_dataset.txt")
pheno_v2<-pheno[,c('X.y','Y','County','trad','barcode_id','Individual.sbarcodenumber.x','Uniqueindividualnumberforthisdate.x')]

both2=merge(pheno_v2,both,by.x='Individual.sbarcodenumber.x',by.y='Barcode')

# climate data
data2=read.delim('~/Dropbox/Amanda_files/ayroles_lab/Turkana_project/phenotypic_data/24Jan19_birth_location/Kenya_precip_2.5min.txt')
data3=read.delim('~/Dropbox/Amanda_files/ayroles_lab/Turkana_project/phenotypic_data/24Jan19_birth_location/Kenya_tavg_2.5min.txt')
data4=read.delim('~/Dropbox/Amanda_files/ayroles_lab/Turkana_project/phenotypic_data/24Jan19_birth_location/Kenya_tmax_2.5min.txt')

data4$avg<-apply(data4[,1:12],1,mean)
data2$avg<-apply(data2[,1:12],1,mean)

# avg in rectangle
data4_t<-subset(data4,x>35 & x<36 & y<4.5 & y>2.5)
data4_t2<-melt(data4_t[,1:14], id.vars=c("x", "y"))

ggplot(data4_t2, aes(x=variable, y=(value* 9/5) + 32)) +geom_boxplot(width=0.5,outlier.shape = NA)+ylim(80,100)+theme_bw(13)+ylab('')
  
data2_t<-subset(data2,x>35 & x<36 & y<4.5 & y>2.5)
data2_t2<-melt(data2_t[,1:14], id.vars=c("x", "y"))

ggplot(data2_t2, aes(x=variable, y=value/10)) + geom_boxplot(width=0.5,outlier.shape = NA)+theme_bw(13)+ylim(0,10)+ylab('')

# plot maps
register_google('AIzaSyBqRWaHZTtsvaSOIwDdHjmTVPD-cXvaKrE')

s <- "element:geometry%7Ccolor:0xf5f5f5&style=element:labels%7Cvisibility:off&style=element:labels.icon%7Cvisibility:off&style=element:labels.text.fill%7Ccolor:0x616161&style=element:labels.text.stroke%7Ccolor:0xf5f5f5&style=feature:administrative%7Celement:geometry%7Cvisibility:off&style=feature:administrative.country%7Celement:geometry.stroke%7Ccolor:0x000000%7Cvisibility:on&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.land_parcel%7Celement:labels.text.fill%7Ccolor:0xbdbdbd&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:poi%7Cvisibility:off&style=feature:poi%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:poi%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:poi.park%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:poi.park%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:road%7Cvisibility:off&style=feature:road%7Celement:geometry%7Ccolor:0xffffff&style=feature:road%7Celement:labels.icon%7Cvisibility:off&style=feature:road.arterial%7Celement:labels.text.fill%7Ccolor:0x757575&style=feature:road.highway%7Celement:geometry%7Ccolor:0xdadada&style=feature:road.highway%7Celement:labels.text.fill%7Ccolor:0x616161&style=feature:road.local%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&style=feature:transit%7Cvisibility:off&style=feature:transit.line%7Celement:geometry%7Ccolor:0xe5e5e5&style=feature:transit.station%7Celement:geometry%7Ccolor:0xeeeeee&style=feature:water%7Celement:geometry%7Ccolor:0xc9c9c9&style=feature:water%7Celement:labels.text.fill%7Ccolor:0x9e9e9e&size=480x360"
map <- get_googlemap(center=c(36,2.25), zoom = 7, scale = 1, style = s)
m<- ggmap(map)

# all samples
data4_v2<-subset(data4,x> 34 &x<38 &y>-0.5 & y<5)
wgs<-subset(both2,type=='WGS')
wgs2<-as.data.frame(table(wgs$X.y,wgs$Y))
wgs2<-subset(wgs2,Freq>0)
names(wgs2)[3]<-'N'
write.table(wgs2,'temp1.txt',sep='\t')
wgs2=read.delim('temp1.txt')

wgs<-subset(both2,type!='WGS')
wgs2<-as.data.frame(table(wgs$X.y,wgs$Y))
wgs2<-subset(wgs2,Freq>0)
names(wgs2)[3]<-'N'
write.table(wgs2,'temp2.txt',sep='\t')
array2=read.delim('temp2.txt')
wgs2=read.delim('temp1.txt')

array2$Sample<-'Array'
wgs2$Sample<-'WGS'
both<-as.data.frame(rbind(array2,wgs2))
m + scale_fill_viridis_c()+geom_jitter(data=both, aes(x = Var1, y = Var2,size=N,color=Sample),width = 0.05, height = 0.05,alpha=0.5)+xlab('Longitude')+ylab('Latitude')+scale_color_manual(values=c('red','steelblue'))+theme_bw(13)+xlim(34,37.5)

# array only with climate
data4_v2<-subset(data4,x> 34 &x<38 &y>-0.5 & y<5)
wgs<-subset(both2,type!='WGS')
m + geom_tile(data=data4_v2, aes(x = x, y = y, fill = avg), alpha=0.75)+scale_fill_viridis_c()+geom_jitter(data=wgs, aes(x = X.y, y = Y),width = 0.05, height = 0.05,alpha=0.1)+xlab('Longitude')+ylab('Latitude')

##########
# water
##########

pheno_past<-subset(pheno,Tribe=='Turkana' & trad=='1_trad')

table(pheno_past$If.you.had.more.water.would.you.drink.more.)
table(pheno_past$How.much.of.your.day.is.spent.finding..carrying..or.retrieving.water.)
table(pheno_past$How.often.do.you.eat.meat.)

y<-c(47,136,152,11)
both<-as.data.frame(y)
both$x<- factor(x, levels=c('Everyday','>2x per week','1-2x per week', 'Rarely'))

ggplot(both, aes(y=y, x=x)) + geom_bar(stat="identity")+theme_bw(13)+coord_flip()+ylab('Frequency')

pheno_past$Specific.gravity[which(pheno_past$Specific.gravity==1.050)]<-1.03
tmp<-as.data.frame(table(pheno_past$Specific.gravity)) 
write.table(tmp,'temp.txt',sep='\t')
tmp<-read.delim('temp.txt')
tmp$Category<-'Hydrated'
tmp$Category[tmp$Var1>1.009]<-'Minimal'
tmp$Category[tmp$Var1>1.024]<-'Significant'

library(RColorBrewer)
cols<-c('red','firebrick','goldenrod1')
ggplot(tmp, aes(y=Freq, x=Var1,fill=Category)) + geom_bar(stat="identity",alpha=0.75)+theme_bw(13)+xlab('Urine specific gravity')+scale_fill_manual(values=cols)+ theme(legend.position="bottom")
    
##########
# diet maps
##########

pheno$milk2<-pheno$milk
pheno$milk2[which(pheno$milk==0)]<-1
pheno$milk2[which(pheno$milk==4)]<-3

tmp<-as.data.frame(table(factor(pheno$milk2),(pheno$Standardized_name)))

library(reshape2)
tmp2 <- as.data.frame(dcast(tmp, Var2 ~ Var1, value.var="Freq"))
write.table(tmp2,'temp.txt',sep='\t')
tmp2<-read.delim('temp.txt')

tmp2<-unique(merge(tmp2,pheno[,c('Standardized_name','X.y','Y')],by.x='Var2',by.y='Standardized_name'))
tmp2$total<-apply(tmp2[,c(2:4)],1,sum)
tmp2<-subset(tmp2,total>10)

names(tmp2)[2:4]<-c('Never_Rarely','Sometimes','Frequently_Always')
library(scatterpie)
m + geom_scatterpie(aes(x=X.y, y=Y), data=tmp2, alpha=0.75,cols=c( "Never_Rarely", "Sometimes",'Frequently_Always'),pie_scale=2.5) + xlim(33.5,38)+scale_fill_brewer(palette='Set1')+theme_bw(13)

##

tmp<-as.data.frame(table(factor(pheno$usg2),(pheno$County),pheno$Wet_season))

ggplot(tmp, aes(fill=Var1, y=Freq, x=Var2)) + geom_bar(position="fill", stat="identity")+facet_wrap(~Var3)