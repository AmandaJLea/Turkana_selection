############
# normalization
############

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/PBMCs/20Aug21_NovaSeq")
library(data.table)
library(limma)
library(edgeR)
library(sva)
library(EMMREML)
library(corrplot)
library(ggplot2)

# raw count data, filter for expressed protein coding genes
data=fread('24Aug21_feature_counts.txt')
protein_coding=read.delim('~/Dropbox/Amanda_files/ayroles_lab/cell_stimulation_exp/Feb2019_NovaSeq/data/protein_coding_hg38.txt')
data2=as.data.frame(subset(data,Gene %in% protein_coding$Gene.stable.ID))

tot_counts<-apply(data2[,-1],2,sum)
remove<-which(tot_counts<250000)
data3<-data2[,-c(remove+1)]
tot_counts<-apply(data3[,-1],2,sum)
norm1=10^6 *(data3[,-1]/tot_counts)
med_TPM<-apply(norm1,1,median)
data4<-data3[which(med_TPM>1),]

# from Audrey - RNASeq metadata
meta=read.delim('RNASeqIndividDataNoCellCount.txt')
count_names=as.data.frame(names(data4[,-1]))
names(count_names)<-'Sample'
count_names$order<-1:dim(count_names)[1]
count_names2=merge(count_names,meta,by='Sample',all.x=T)

keep<-subset(count_names2,  Environment!='NA' & Age!='NA' & Gender!='NA' )
design=model.matrix(~ Gender + Age + Environment, data = keep)

#normalize with voom
 dge <- DGEList(counts=data4[,c(keep$Sample)])
 dge <- calcNormFactors(dge)
 v <- voom(dge,design, plot = TRUE)

# sent to Kristina
# write.table(v$E,'1Jul22_normalized_expression.txt',row.names=F,sep='\t',quote=F)

# control for sequencing flow cell batch effects
v_combat = ComBat(dat=as.matrix(v$E), batch=keep$Lane, mod=design, par.prior=TRUE)
v$E=v_combat
rownames(v$E)<-data4$Gene

############
# add cell count and genetic data
############

# cell counts
cells=read.csv('~/Dropbox/Amanda_files/ayroles_lab/Turkana_project/early_life/data/Data Entry_Turkana_latest_2020 - Cell count.csv')
cells$barcode_id<-paste(cells$Unique.barcode,cells$Individual.number.for.this.date,sep='_')
cells$count_tot<-apply(cells[,c(6:10)],1,sum)

cells$neut<-cells$Segmented.neutrophils/cells$count_tot
cells$mono<-cells$Monocytes/cells$count_tot
cells$lymp<-cells$Lymphocytes/cells$count_tot
cells$eos<-cells$Eosinophils/cells$count_tot
cells$baso<-cells$Basophils/cells$count_tot
cells<-subset(cells,count_tot>99)

pca_cells<-prcomp(cells[,c('neut','mono','lymp','eos','baso')])
pr.var <- pca_cells$sdev^2
pve <- pr.var / sum(pr.var)
cells$cell_PC1<-pca_cells$x[,1]
cells2<-as.data.frame(aggregate(cells$cell_PC1 ~ cells$barcode_id,FUN=mean))
names(cells2)<-c('barcode_id','cell_PC1')
cells3<-unique(merge(cells2,cells[,c('Unique.barcode','Individual.number.for.this.date','barcode_id')],by='barcode_id'))

# genetic data
eigenvec=read.delim('/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/PBMCs/20Aug21_NovaSeq/24Jun22_analyses/merge_8batches_filt.eigenvec.txt',header=T)
eigenvec$geno_matrix_id<-1:dim(eigenvec)[1]
rel=read.delim('/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/PBMCs/20Aug21_NovaSeq/24Jun22_analyses/merge_8batches_filt.rel',header=F)

# merge
keep$rna_matrix_id<-1:dim(keep)[1]
keep2<-merge(keep,eigenvec[,-c(1,2,10:14)],all.x=T,by.x='SampleNumber',by.y='Barcode')
keep3<-merge(keep2,cells3,all.x=T,by.x='SampleNumber',by.y='Unique.barcode')

# with cell type
keep4_wcell<-subset(keep3,Environment!='NA' & Age!='NA' & Gender!='NA' &X.2> -100 & cell_PC1>-100)
dup<-subset(as.data.frame(table(keep4_wcell$SampleNumber)),Freq>1)

keep4_wcell$cell_correctID<-0
keep4_wcell$cell_correctID[-which(keep4_wcell$SampleNumber %in% dup$Var1)]<-1
tmp<-subset(keep4_wcell,SampleNumber %in% dup$Var1)
tmp<-subset(tmp,Individual.number.for.this.date.x==Individual.number.for.this.date.y)
keep4_wcell$cell_correctID[which(keep4_wcell$barcode_id %in% tmp$barcode_id)]<-1
keep4_wcell<-subset(keep4_wcell,cell_correctID==1)

library(ggplot2)
ggplot(keep4_wcell, aes(x=X.2, y=X.3, color=Lifestyle)) +geom_point()+ylim(-0.05,0.05)

############
# PCA
############

pca_cells<-prcomp(t(v$E[,keep4_wcell$rna_matrix_id]))
pr.var <- pca_cells$sdev^2
pve <- pr.var / sum(pr.var)
keep4_wcell$exp_PC1<-pca_cells$x[,1]
keep4_wcell$exp_PC2<-pca_cells$x[,2]

ggplot(keep4_wcell, aes(x=exp_PC1, y=exp_PC2, color=Environment)) +geom_point()
ggplot(keep4_wcell, aes(x=exp_PC1, y=exp_PC2, color=Lifestyle)) +geom_point()

cut<-median(pca_cells$x[,1])+sd(pca_cells$x[,1])*4

# outliers removed and normalize samp size
keep4_wcell<-subset(keep4_wcell,exp_PC1<cut )

keep4_wcell$keep<-0
keep4_wcell$keep[which(keep4_wcell$Lifestyle!='Rural')]<-1
x<-length(keep4_wcell$keep[which(keep4_wcell$Lifestyle=='Urban')])
keep4_wcell$keep[which(keep4_wcell$Lifestyle=='Urban')]<-sample(c(rep(0,x-115),rep(1,115)))

keep4_wcell<-subset(keep4_wcell,keep==1)

############
# model with genetic data and with cell type
############

to_use_metadata=keep4_wcell
to_use_data=v$E[,keep4_wcell$rna_matrix_id]

design=model.matrix(~  Age + Environment +cell_PC1 + X.2, data = to_use_metadata)
k_matrix<-rel[to_use_metadata$geno_matrix_id,to_use_metadata$geno_matrix_id]

out<-matrix(nrow=dim(to_use_data)[1],ncol=5)
out_beta<-matrix(nrow=dim(to_use_data)[1],ncol=5)

for (i in 1:dim(to_use_data)[1]) {
emma=emmreml(y=(to_use_data[i,]),X=design,Z=as.matrix(diag(dim(k_matrix)[2])),K=as.matrix(k_matrix),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
	out[i,]<-emma$pvalbeta[,8]
	out_beta[i,]<-emma$betahat
print(i) }

par(mfrow=c(2,3))
names<-c('Age','Env','Celltype_PC','Geno_PC1')
for (i in 2:5){
	hist(out[,i],main=names[i-1])
}

out2<-as.data.frame(out)
out2$env_qvalue<-qvalue(out[,3])$qvalues
out2$env_beta<-out_beta[,3]
out2$geneID<-rownames(v$E)
write.table(out2,'20Sep22_LMM_pvals.txt',row.names=F,sep='\t')
write.table(out_beta,'20Sep22_LMM_betas.txt',row.names=F,sep='\t')
write.table(to_use_metadata,'20Sep22_LMM_metadata.txt',row.names=F,sep='\t')

##

ens<-c('ENSG00000204388','ENSG00000204389','ENSG00000173110','ENSG00000132002','ENSG00000169245','ENSG00000080824','ENSG00000151929','ENSG00000144381','ENSG00000140403','ENSG00000086061','ENSG00000162783','ENSG00000096384','ENSG00000142168')
name<-c('HSPA1B (Hsp70)','HSPA1A (Hsp70)','HSPA6 (Hsp70)','DNAJB1 (Hsp40)','CXCL10','HSP90AA1 (Hsp90)','BAG3','HSPD1 (Hsp60)','DNAJA4 (Hsp40)','DNAJA1 (Hsp40)','IER5','HSP90AB1 (Hsp90)','SOD1')

to_use_metadata=read.delim('20Sep22_LMM_metadata.txt')
v_keep<-as.data.frame(v$E[ens,to_use_metadata$rna_matrix_id])
v_keep$gene<-name
library(reshape2)
v_keep2<-melt(v_keep,id.vars=c('gene'))
v_keep2<-merge(v_keep2,to_use_metadata[,c('Environment','Gender','Sample')],by.x='variable',by.y='Sample')
v_keep2<-subset(v_keep2,gene %in% c('HSPA1B (Hsp70)','HSPA1A (Hsp70)','HSPA6 (Hsp70)','DNAJB1 (Hsp40)','HSP90AA1 (Hsp90)','HSPD1 (Hsp60)','DNAJA4 (Hsp40)','DNAJA1 (Hsp40)','HSP90AB1 (Hsp90)'))

data_summary <- function(x) {
   m <- mean(x)
   ymin <- m-sd(x)
   ymax <- m+sd(x)
   return(c(y=m,ymin=ymin,ymax=ymax))
}

library(ggplot2)
ggplot( aes(x=reorder(gene, value), y=value, fill=Environment),data=v_keep2) +
    geom_violin(position=position_dodge(1),width=1.5) +coord_flip()+theme_bw(13)+scale_fill_manual(values=c('steelblue','goldenrod'))+ stat_summary(fun.data=data_summary,position=position_dodge(1))
    
to_use_metadata=read.delim('20Sep22_LMM_metadata.txt')
v_keep<-as.data.frame(v$E[ens,to_use_metadata$rna_matrix_id])
v_keep$gene<-name
library(reshape2)
v_keep2<-melt(v_keep,id.vars=c('gene'))
v_keep2<-merge(v_keep2,to_use_metadata[,c('Lifestyle','Gender','Sample')],by.x='variable',by.y='Sample')
v_keep2<-subset(v_keep2,gene %in% c('HSPA1B (Hsp70)','HSPA1A (Hsp70)','HSPA6 (Hsp70)','DNAJB1 (Hsp40)','HSP90AA1 (Hsp90)','HSPD1 (Hsp60)','DNAJA4 (Hsp40)','DNAJA1 (Hsp40)','HSP90AB1 (Hsp90)'))

ggplot( aes(x=reorder(gene, value), y=value, fill=Lifestyle),data=v_keep2) +
    geom_violin(position=position_dodge(1),width=1.5) +coord_flip()+theme_bw(15)+scale_fill_manual(palette='Set1')+ stat_summary(fun.data=data_summary,position=position_dodge(1))

ggplot( aes(x=Lifestyle, y=value, fill=Lifestyle),data=v_keep2) +
    geom_violin(position=position_dodge(1)) +theme_bw(13)+scale_fill_brewer(palette='Set1')+ stat_summary(fun.data=data_summary,position=position_dodge(1))+facet_wrap(~gene,scales='free_y')+ylab('Gene expression level')+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

############
# GSEA
############

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/PBMCs/20Aug21_NovaSeq")
out_coef=read.delim('20Sep22_LMM_pvals.txt')

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# using abs value
original_gene_list <- -1*out_coef$env_beta
names(original_gene_list) <- out_coef$geneID
gene_list<-na.omit(original_gene_list)
gene_list = sort( (gene_list), decreasing = TRUE)

gse_res <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             minGSSize = 10, maxGSSize = 500,
             pvalueCutoff = 1, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

gse_res2<-as.data.frame(gse_res)
write.table(gse_res2[,1:8],('20Sep22_GSEA_DE_genes.txt'),sep='\t',row.names=F)

x2 <- pairwise_termsim(gse_res)
emapplot(x2)

############
# plot
############

library(GOplot)
circ <- circle_dat(gse_res$david, EC$core_enrichment)

