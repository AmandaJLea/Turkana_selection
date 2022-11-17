#############
# RUN IN BASH
#############

module load bcftools

# extract pops of interest
bcftools view -O z -o /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.YRI_LWK_CEU_MKK.vcf.gz --samples-file YRI_LWK_CEU_MKK.txt /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz --force-samples

bcftools index /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.YRI_LWK_CEU_MKK.vcf.gz

# harmonize
path_harmonizer=~/programs/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar
java -Xmx50g -jar $path_harmonizer --inputType PLINK_BED --input ~/Turkana/genotyping_Jul22/merge_9batches_filt1 --ref /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.YRI_LWK_CEU_MKK.vcf.gz --output ~/Turkana/genotyping_Jul22/merge_9batches_filt1_harm  --outputType PLINK_BED --update-reference-allele --ambiguousSnpFilter --update-id

# merge array data
~/programs/plink_1.90 --vcf /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.YRI_LWK_CEU_MKK.vcf.gz --bmerge ~/Turkana/genotyping_Jul22/merge_9batches_filt1_harm --geno 0.25 --maf 0.01  --indep-pairwise 50 20 0.8 --out ~/Turkana/genotyping_Jul22/merge_9batches_filt1_w1000G --allow-extra-chr --chr 1-22

# exclude
~/programs/plink_1.90 --vcf /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.YRI_LWK_CEU_MKK.vcf.gz --out temp1 --make-bed --exclude /Genomics/grid/users/alea/Turkana/genotyping_Jul22/merge_9batches_filt1_w1000G.missnp --allow-extra-chr --chr 1-22

~/programs/plink_1.90 --bfile ~/Turkana/genotyping_Jul22/merge_9batches_filt1_harm --out temp2 --make-bed --exclude /Genomics/grid/users/alea/Turkana/genotyping_Jul22/merge_9batches_filt1_w1000G.missnp

# merge
~/programs/plink_1.90 -bfile temp1 --bmerge temp2 --geno 0.01 --maf 0.05 --indep-pairwise 50 20 0.8 --out ~/Turkana/genotyping_Jul22/merge_9batches_filt1_w1000G --chr 1-22

# extract
~/programs/plink_1.90 -bfile ~/Turkana/genotyping_Jul22/merge_9batches_filt1_w1000G --out ~/Turkana/genotyping_Jul22/merge_9batches_filt1_w1000G_LDfilt --make-bed --extract /Genomics/grid/users/alea/Turkana/genotyping_Jul22/merge_9batches_filt1_w1000G.prune.in

# run ADMIXTURE
K=3
~/programs/admixture_linux-1.3.0/admixture --cv ~/Turkana/genotyping_Jul22/merge_9batches_filt1_w1000G_LDfilt.bed $K --seed=$RANDOM | tee ~/Turkana/genotyping_Jul22/merge_9batches_filt1_w1000G_LDfilt.$K  

#############
# RUN IN R
#############

MKK=c('NA21295','NA21296','NA21297','NA21298','NA21299','NA21300','NA21301','NA21302','NA21303','NA21304','NA21305','NA21306','NA21307','NA21308','NA21309','NA21310','NA21311','NA21312','NA21313','NA21314','NA21316','NA21317','NA21318','NA21319','NA21320','NA21321','NA21333','NA21334','NA21335','NA21336','NA21337','NA21338','NA21339','NA21340','NA21341','NA21344','NA21351','NA21352','NA21353','NA21355','NA21356','NA21357','NA21358','NA21359','NA21360','NA21361','NA21362','NA21363','NA21364','NA21365','NA21366','NA21367','NA21368','NA21369','NA21370','NA21371','NA21378','NA21379','NA21380','NA21381','NA21382','NA21383','NA21384','NA21385','NA21386','NA21387','NA21388','NA21389','NA21390','NA21391','NA21399','NA21400','NA21401','NA21402','NA21403','NA21404','NA21405','NA21407','NA21408','NA21409','NA21410','NA21414','NA21415','NA21416','NA21417','NA21418','NA21419','NA21420','NA21421','NA21423','NA21424','NA21425','NA21433','NA21434','NA21435','NA21436','NA21438','NA21439','NA21440','NA21441','NA21442','NA21443','NA21444','NA21447','NA21448','NA21449','NA21450','NA21451','NA21452','NA21453','NA21454','NA21455','NA21456','NA21457','NA21458','NA21473','NA21474','NA21475','NA21476','NA21477','NA21478','NA21479','NA21480','NA21485','NA21486','NA21487','NA21488','NA21489','NA21490','NA21491','NA21492','NA21493','NA21494','NA21509','NA21510','NA21511','NA21512','NA21513','NA21514','NA21515','NA21517','NA21518','NA21519','NA21520','NA21521','NA21522','NA21523','NA21524','NA21525','NA21526','NA21527','NA21528','NA21529','NA21530','NA21573','NA21574','NA21575','NA21576','NA21577','NA21578','NA21579','NA21580','NA21581','NA21582','NA21583','NA21584','NA21586','NA21587','NA21588','NA21596','NA21597','NA21599','NA21600','NA21601','NA21608','NA21609','NA21611','NA21613','NA21614','NA21615','NA21616','NA21617','NA21618','NA21619','NA21620','NA21622','NA21631','NA21632','NA21634','NA21635','NA21636','NA21647','NA21648','NA21649','NA21650','NA21678','NA21679','NA21682','NA21683','NA21685','NA21686','NA21689','NA21693','NA21716','NA21717','NA21718','NA21719','NA21722','NA21723','NA21732','NA21733','NA21735','NA21737','NA21738','NA21739','NA21740','NA21741','NA21742','NA21743','NA21744','NA21768','NA21769','NA21770','NA21774','NA21775','NA21776','NA21782','NA21784','NA21785','NA21786','NA21825','NA21826')

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures")

meta1=read.delim("15Jul22_all_genotyped_samples.txt",header=T)
meta2=read.delim('/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/25Aug19_NovaSeq_FINAL/CURRENT_pipeline/data/20130606_g1k.ped.txt')
data=read.table('merge_9batches_filt1_w1000G_LDfilt.3.Q',header=F,sep=' ')
samps=read.delim('merge_9batches_filt1_w1000G_LDfilt.fam',header=F,sep=' ')
data$samp<-paste(samps$V2)

data2<-merge(data,meta1,by.x='samp',by.y='FileID',all.x=T)
data2$Group[which(data2$samp %in% subset(meta2,Population=='CEU')$Alt.Individual.ID)]<-'CEU'
data2$Group[which(data2$samp %in% subset(meta2,Population=='LWK')$Alt.Individual.ID)]<-'LWK'
data2$Group[which(data2$samp %in% subset(meta2,Population=='YRI')$Alt.Individual.ID)]<-'YRI'
data2$Group[which(data2$samp %in% MKK)]<-'MKK'

# downsample large groups
downsample<-c('CEU','Turkana','LWK','YRI','MKK')
data2$keep<-1
data2$keep[which(data2$Group %in% downsample)]<-0

for (f in downsample){
tmp<-subset(data2,Group==f)
data2$keep[which(data2$samp %in% tmp$samp)]<-sample(c(rep(1,20),rep(0,dim(tmp)[1] -20))) }

data3<-subset(data2,keep==1)
data3<-subset(data3,Group!='NA')

library(tidyr)

data_long <- gather(data3[,c('samp','Group','V1','V2','V3')], admixture_comp, admixture_prop, V1:V3, factor_key=TRUE)
data_long<-data_long[order(data_long$samp,data_long$admixture_comp),]
data_long$admixture_comp2<-c(3,2,1)
data_long<-data_long[order(data_long$samp,data_long$admixture_comp2),]

library(RColorBrewer)
cols<-brewer.pal(3,'Paired')

library(ggplot2)
ggplot(data_long, aes(factor(samp), admixture_prop, fill = factor(admixture_comp2))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~Group, scales = "free", space = "free") + labs(x = "Individuals", title = "K=3", y = "Ancestry component") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank() ) + theme(legend.position = "none")+ scale_fill_manual(values=c(cols[1],cols[2],cols[3]))+ theme(strip.text.y = element_text(angle = 180))

