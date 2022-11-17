#########
# enrichment of candidate genes in GWAS hits
#########

cand<-c(2,2,2,2,2,2,2,2,2,2,3,3,4,4)
gwas<-c(382,2132,5631,1041,1558,2964,692,2177,3336,2145,1639,1857,1343,939)

pval<-c()
estimate<-c()

for (i in 1:length(cand)){
mat<-matrix(c(19430-(8-cand[i])-gwas[i], 
gwas[i]-cand[i],
8-cand[i],
cand[i]),ncol=2)

estimate<-c(estimate,fisher.test(mat)$estimate)
pval<-c(pval,fisher.test(mat)$p.value)}

out<-as.data.frame(cbind(pval,estimate,p.adjust(pval,method='BH')))

#########
# enrichment of candidate genes in particular phenotypes
#########

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures")

data=read.delim('9Sep22_trait_associations.txt')
data$id<-1:dim(data)[1]
genes<-unique(data$Gene)
data$top<-0

for (i in 1:length(genes)){
tmp=subset(data,Gene==genes[i] & Trait!='protein' & Trait!='blood protein' & Trait!='self reported educational attainment' & Trait!='mathematical ability' & Trait!='neuroticism' & Trait!='neuroimaging'& Trait!='brain volume' & Trait!='brain')
tmp=tmp[order(tmp$overallAssociationScore,decreasing = TRUE),]
data$top[which( data$id %in% tmp$id[1:5])]<-1 }

library(ggplot2)
tmp<-subset(data,top==1 &  Trait!='protein' & Trait!='blood protein' & Trait!='self reported educational attainment' & Trait!='mathematical ability' & Trait!='neuroticism' & overallAssociationScore>0.25 & Trait!='neuroimaging' & Trait!='brain volume' & Trait!='brain')
tmp$stc1<-0
tmp$stc1[which(tmp$Gene=='STC1')]<-1

library(RColorBrewer)
cols<-brewer.pal(8,'Blues')
cols2<-c(cols[2:6],'red',cols[7:8])

ggplot(data=tmp, aes(x=reorder(Trait, overallAssociationScore), y=overallAssociationScore,color=Gene,shape=factor(stc1))) +geom_linerange(aes(ymin = 0, ymax = overallAssociationScore),size=1.1)+geom_point(stat="identity",size=3)+coord_flip()+theme_bw(12)+ylab('Association Score')+xlab('')+scale_color_manual(values=cols2)+scale_shape_manual(values=c(20,8))

