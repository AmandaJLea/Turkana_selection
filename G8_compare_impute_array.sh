##########
# RUN IN BASH
##########

cd /Genomics/grid/users/alea/Turkana/genotyping_Jul22

# export chr1 - imputed WGS
~/programs/plink_1.90 --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.ALL.hg19_GQ10_plink --chr 1 --recode A-transpose --out impute1panel.1.ALL.hg19_GQ10_plink_chr1 --hwe 0.000001

# export chr1 - array
~/programs/plink_1.90 --bfile ~/Turkana/genotyping_Jul22/merge_9batches_filt1 --remove kingunrelated_toberemoved3.txt --recode A-transpose --out merge_9batches_filt1_chr1 

~/programs/plink_1.90 --bfile ~/Turkana/genotyping_Jul22/merge_9batches_filt1 --remove kingunrelated_toberemoved3.txt --make-just-fam --out merge_9batches_filt1_chr1 


##########
# RUN IN R
##########

library(data.table)
wgs=fread('impute1panel.1.ALL.hg19_GQ10_plink_chr1.traw')
array=fread('merge_9batches_filt1_chr1.traw')

wgs_names<-read.delim('/scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.ALL.hg19_GQ10_plink.fam',header=F,sep=' ')
array_names<-read.delim('merge_9batches_filt1_chr1.fam',header=F,sep=' ')

wgs_names$id<-1:dim(wgs_names)[1]
array_names$id<-1:dim(array_names)[1]

both<-merge(wgs_names,array_names,by.x='V1',by.y='V2')
both<-subset(both,V1<600 | V1>1100)

# included in final WGS dataset
keep<-c(1090,1091,1095,1098,1101,1103,1105,1106,1212,1242,1275,1310,1311,1170,1182,1185,1220,1224,1244,1267,1274,596,598,603,608,609,613,614,1196,1263,1266,1269,1282,1297,1299,1312,1035,1036,1038,1041,1044,1046,1047,1048,1051,1052,1055,1056,1058,1060,1063,1070,971,974,980,981,985,987,997,16,22,28,32,49,62,64,160,161,164,166,168,170,172,173,174,175,176,177,179,184,188,189,191,193,197,207,208,210,214,215,216,217,220,221,225,227,229,230,231,236,240,241,244,248,250,251,254,255,256,257,264,267,269,270,276,285,286,287,290,291,299,307,318,319,324,337,340,351,353,358,359,372,375,376,378,385,386,388,389,390,399,403,405,419,436,442,446,454,471,481,485,486,490,492,493,496,497,498,499,500,504,589,593,616,624,636,645,648,652,654,659,663,664,665,672,677,681,682,693,694,696,699,703,705,706,753,754,760,775,791,799,801,803,805,809,817,820,821,822,823,824,832,833,838,849,850,852,854,856,858,859,863,864,876,879,885,886,887,890,895,898,899,900,901,902,903,904,905,906,908,910,911,914,918,922,923,929,931,932,938,945,950,953,957,966,968,970,981,992,1003,1004,1007,1008,1012,1017,1025,1039,1042,1044,1057,1177)
both<-subset(both,V1 %in% keep)

cor<-c()

for (i in 1:dim(both)[1]){
tmp1<-wgs[,c(1:6,both$id.x[i]+6),with=F]
tmp2<-array[,c(1:6,both$id.y[i]+6),with=F]
tmp3<-merge(tmp1,tmp2,by='POS')
tmp4<-as.data.frame(subset(tmp3,COUNTED.x==COUNTED.y & ALT.x==ALT.y))
cor<-c(cor, cor.test(tmp4[,7],tmp4[,13],method='spearman')$estimate) }

