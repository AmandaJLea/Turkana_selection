##########
# RUN IN BASH
###########

# phenotype manifest, from Pan UKBB: https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=30994804

##########
# consider a window as associated with a trait if it included a SNP with a genome-wide significant association with this trait (p<10^-8)
###########

# cd /Genomics/ayroleslab2/alea/ref_genomes/public_datasets/UKBB_pan/
ls *bgz > files.txt

library(data.table)
library(stringr)

files=read.delim('files.txt',header=F)
files2=str_split_fixed(files$V1, "-", 4)

for (i in 1:38){
x=paste("zcat ",files$V1[i],sep="")
y=fread(x)
sig=subset(y,pval_meta<10^-8 & low_confidence_MID=='FALSE')

sig$loc2<-sig$pos+1
write.table(sig[,c('chr','pos','loc2','pval_meta')],paste('trait-',files2[i,2],'-sig.bed',sep=''),row.names=F,sep='\t',col.names=F,quote=F)
}

for (i in 39:40){
x=paste("zcat ",files$V1[i],sep="")
y=fread(x)
sig=subset(y,pval_meta<10^-8 & low_confidence_AFR=='FALSE')
sig$loc2<-sig$pos+1
write.table(sig[,c('chr','pos','loc2','pval_meta')],paste('trait-',files2[i,2],'-sig.bed',sep=''),row.names=F,sep='\t',col.names=F,quote=F)
}

#

windows=/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_v0_Homo_sapiens_assembly19.50kb_windows.nooverlap_sort.txt
module load bedtools

for f in trait*sig.bed; do bedtools intersect -a $windows -b $f -u > windows.$f; done

cat windows.trait*sig.bed | awk '{print $1,'_',$2}' | sort | uniq -c | wc -l #18255