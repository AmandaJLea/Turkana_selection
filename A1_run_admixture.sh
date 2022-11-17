#############
# RUN IN BASH
#############

# prep and merge
plink=/Genomics/grid/users/alea/programs/plink_1.90
path_resources=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets
samples=YRI_LWK_CEU.txt

# make bed for 1000G
$plink --vcf $path_resources/shapeit_files/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out $path_resources/shapeit_files/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.YRI_LWK_CEU --indep-pairwise 50 20 0.8 --maf 0.01 --keep YRI_LWK_CEU.txt

# change rsID
awk '{OFS="\t";print $1,$1":"$4,$3,$4,$5,$6}' $path_resources/shapeit_files/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.YRI_LWK_CEU.bim > temp.bim
mv temp.bim $path_resources/shapeit_files/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.YRI_LWK_CEU.bim

# first pass merge
$plink --bfile /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/impute1panel.1.merged_LDfilt_plink --bmerge $path_resources/shapeit_files/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.YRI_LWK_CEU --out temp1

# filter 1000G
$plink --bfile $path_resources/shapeit_files/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.YRI_LWK_CEU --make-bed --exclude temp1.missnp --out $path_resources/shapeit_files/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.YRI_LWK_CEU_v2

# filter Turkana
$plink --bfile $LDfilt --make-bed --exclude temp1.missnp --out temp2

# final merge
$plink --bfile $path_resources/shapeit_files/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.YRI_LWK_CEU_v2 --make-bed --out temp3 --bmerge temp2 --geno 0.25 --maf 0.01 --indep-pairwise 50 20 0.8

$plink --bfile temp3 --make-bed --out impute1panel.1.merged_LDfilt_forADMIX --exclude temp3.prune.out 

# run ADMIXTURE
for K in {2..7} ; do ~/programs/admixture_linux-1.3.0/admixture --cv impute1panel.1.merged_LDfilt_forADMIX.bed $K --seed=$RANDOM | tee impute1panel.1.merged_LDfilt_forADMIX.${K}.$RANDOM.out  ; done

#############
# RUN IN R
#############

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures")

cv=read.delim('29Aug22_admixture_CV.txt')

library(ggplot2)
  ggplot( aes(x=factor(K), y=CV),data=cv) +  geom_boxplot(outlier.shape = NA) +  geom_jitter(color="black", alpha=0.4) + theme_bw(13)+xlab('K')
######### K = 3

meta1=read.delim("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data/impute1panel.1.merged_LDfilt_plink.fam_winfo.txt",header=F)
meta1$samp<-paste(meta1$V1,meta1$V2,sep='_')
meta2=read.delim('/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data/20130606_g1k.ped.txt')
data=read.table('impute1panel.1.merged_LDfilt_forADMIX.3.Q',header=F,sep=' ')
samps=read.delim('impute1panel.1.merged_LDfilt_forADMIX.fam',header=F,sep=' ')
data$samp<-paste(samps$V1,samps$V2,sep='_')
data$samp2<-paste(samps$V1)

data2<-merge(data,meta1,by='samp',all.x=T)
data2$V3.y[which(data2$samp2 %in% subset(meta2,Population=='CEU')$Alt.Individual.ID)]<-'CEU'
data2$V3.y[which(data2$samp2 %in% subset(meta2,Population=='LWK')$Alt.Individual.ID)]<-'LWK'
data2$V3.y[which(data2$samp2 %in% subset(meta2,Population=='YRI')$Alt.Individual.ID)]<-'YRI'

# downsample large groups
downsample<-c('CEU','Turkana','LWK','YRI')
data2$keep<-1
data2$keep[which(data2$V3.y %in% downsample)]<-0

for (f in downsample){
tmp<-subset(data2,V3.y==f)
data2$keep[which(data2$samp %in% tmp$samp)]<-sample(c(rep(1,10),rep(0,dim(tmp)[1] -10))) }

data3<-subset(data2,keep==1)

library(tidyr)

data_long <- gather(data3[,c('samp','V3.y','V1.x','V2.x','V3.x')], admixture_comp, admixture_prop, V1.x:V3.x, factor_key=TRUE)
data_long<-data_long[order(data_long$samp),]

library(ggplot2)
ggplot(data_long, aes(factor(samp), admixture_prop, fill = factor(admixture_comp))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~V3.y, scales = "free", space = "free") + labs(x = "Individuals", title = "K=3", y = "Ancestry component") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank() ) + theme(legend.position = "none")+ scale_fill_brewer(palette="Paired")+ theme(strip.text.y = element_text(angle = 180))

# check meta data
tmp<-subset(data2,V3.y=='Turkana')
summary(lm(tmp$V1.x~tmp$V5+tmp$V4))
summary(lm(tmp$V2.x~tmp$V5+tmp$V4))
summary(lm(tmp$V3.x~tmp$V5+tmp$V4))

boxplot(tmp$V2.x~tmp$V5*tmp$V4)

######### K =4

meta1=read.delim("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data/impute1panel.1.merged_LDfilt_plink.fam_winfo.txt",header=F)
meta1$samp<-paste(meta1$V1,meta1$V2,sep='_')
meta2=read.delim('/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data/20130606_g1k.ped.txt')
data=read.table('impute1panel.1.merged_LDfilt_forADMIX.4.Q',header=F,sep=' ')
samps=read.delim('impute1panel.1.merged_LDfilt_forADMIX.fam',header=F,sep=' ')
data$samp<-paste(samps$V1,samps$V2,sep='_')
data$samp2<-paste(samps$V1)

data2<-merge(data,meta1,by='samp',all.x=T)
data2$V3.y[which(data2$samp2 %in% subset(meta2,Population=='CEU')$Alt.Individual.ID)]<-'CEU'
data2$V3.y[which(data2$samp2 %in% subset(meta2,Population=='LWK')$Alt.Individual.ID)]<-'LWK'
data2$V3.y[which(data2$samp2 %in% subset(meta2,Population=='YRI')$Alt.Individual.ID)]<-'YRI'

# downsample large groups
downsample<-c('CEU','Turkana','LWK','YRI')
data2$keep<-1
data2$keep[which(data2$V3.y %in% downsample)]<-0

for (f in downsample){
tmp<-subset(data2,V3.y==f)
data2$keep[which(data2$samp %in% tmp$samp)]<-sample(c(rep(1,10),rep(0,dim(tmp)[1] -10))) }

data3<-subset(data2,keep==1)

library(tidyr)
data_long <- gather(data3[,c('samp','V3.y','V1.x','V2.x','V3.x','V4.x')], admixture_comp, admixture_prop, V1.x:V4.x, factor_key=TRUE)
data_long<-data_long[order(data_long$samp),]

library(ggplot2)
ggplot(data_long, aes(factor(samp), admixture_prop, fill = factor(admixture_comp))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~V3.y, scales = "free", space = "free") + labs(x = "Individuals", title = "K=4", y = "Ancestry component") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank() ) + theme(legend.position = "none")+ scale_fill_brewer(palette="Paired")+ theme(strip.text.y = element_text(angle = 180))
    

