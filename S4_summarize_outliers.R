########
# RUN IN R
########

# BayEnv results
i=1
data1=read.delim(paste('8Jun20_BayEnv_results_chr',i,'.txt',sep=''))
data2=read.delim(paste('8Jun20_BayEnv_region_results_chr',i,'.txt',sep=''))

for (i in 2:22){
tmp1=read.delim(paste('8Jun20_BayEnv_results_chr',i,'.txt',sep=''))
tmp2=read.delim(paste('8Jun20_BayEnv_region_results_chr',i,'.txt',sep=''))
names(tmp1)<-names(data1)
names(tmp2)<-names(data2)
data1=rbind(data1,tmp1)
data2=rbind(data2,tmp2)
print(i)
}

# PBS and iHS results summarized by region
both2_PBS=read.delim('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/PBS_files/21May20_PBS_ALLchr_outlier_count.txt',sep=' ')
both2_iHS=read.delim('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/rehh_files/26May20_iHS_ALLchr_outlier_count.txt',sep=' ')

all=merge(both2_PBS,both2_iHS,by='region',all.x=T)
all=merge(all,data2,by='region',all.x=T)

genes=read.delim('/Genomics/ayroleslab2/alea/turkana_wgs/joint_vcfs/rehh_files/26May20_iHS_test_gene_info.bed',header=F)
genes$region=paste(genes$V1,genes$V2,genes$V3,sep='_')
both2=merge(all,genes,by='region')

# iHS outlier region cutoff
quantile(all$outlier_SNPs.y,seq(0,1,0.01),na.rm=T)[100]

# PBS outlier region cutoff
quantile(all$outlier_SNPs.x,seq(0,1,0.01))[100]

# BayEnv outlier region cutoff
quantile(all$n_outliers,seq(0,1,0.01))[100]

write.table(subset(both2,n_outliers>5 & outlier_SNPs.x>5 & outlier_SNPs.y>14),'outlier_regions.txt',row.names=F,sep='\t')
write.table(all,'outliers_by_regions.txt',row.names=F,sep='\t')

# manhattan plots

setwd("/Users/alea/Dropbox/Amanda_files/ayroles_lab/Turkana_project/selection_manuscript/For_figures")
data=read.delim('outliers_by_regions.txt')
data$dummy<-1:dim(data)[1]
outliers<-read.delim('outlier_regions.txt')
data$region<-paste(data$chr,data$loc1,data$loc2,sep='_')
highlight<-data$dummy[which(data$region %in% outliers$region)]

library(RColorBrewer)
col<-brewer.pal(12,'Paired')
pdf(file='manhattan_plots.pdf')
par(mfrow=c(3,1))
# PBS
manhattan_mod2(data,chr='chr',bp='loc1',p='outlier_SNPs.x',snp='dummy',logp = FALSE,suggestiveline = 6, genomewideline = 6,ylim=c(0,22),highlight=highlight,ylab='PBS',col = c("gray80",  "gray50"))

# BayEnv
manhattan_mod2(data,chr='chr',bp='loc1',p='n_outliers',snp='dummy',logp = FALSE,suggestiveline = 6, genomewideline = 6,ylim=c(0,22),highlight=highlight,ylab='XtX',col = c("gray80",  "gray50"))

# outlier_SNPs.y
manhattan_mod2(data,chr='chr',bp='loc1',p='outlier_SNPs.y',snp='dummy',logp = FALSE,suggestiveline = 15, genomewideline = 15,ylim=c(0,110),highlight=highlight,ylab='iHS',col = c("gray80",  "gray50"))
dev.off()

# manhattan function

manhattan_mod2<-function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
    "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
    genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
    annotatePval = NULL, annotateTop = TRUE, ...) 
{
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x))) 
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) 
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) 
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x))) 
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]])) 
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) 
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) 
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    if (!is.null(x[[snp]])) 
        d = transform(d, SNP = x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
        d$logp <- -log10(d$P)
    }
    else {
        d$logp <- d$P
    }
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        d$pos = d$BP
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position")
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index == 
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
                  lastbase
            }
            ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                i, ]$pos))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
            ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
        names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            }
            else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        }
        else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logp, col = col[icol], pch = 20, ...))
            icol = icol + 1
        }
    }
    if (suggestiveline) 
        abline(h = suggestiveline, col = 'red',lwd=2,lty=2)
  #  if (genomewideline) 
  #      abline(h = genomewideline, col = 'red',lwd=2)
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = 'red', cex=2,pch = 20, 
            ...))
    }
    if (!is.null(annotatePval)) {
        topHits = subset(d, P <= annotatePval)
        par(xpd = TRUE)
        if (annotateTop == FALSE) {
            with(subset(d, P <= annotatePval), textxy(pos, -log10(P), 
                offset = 0.625, labs = topHits$SNP, cex = 0.45), 
                ...)
        }
        else {
            topHits <- topHits[order(topHits$P), ]
            topSNPs <- NULL
            for (i in unique(topHits$CHR)) {
                chrSNPs <- topHits[topHits$CHR == i, ]
                topSNPs <- rbind(topSNPs, chrSNPs[1, ])
            }
            textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
                labs = topSNPs$SNP, cex = 0.5, ...)
        }
    }
    par(xpd = FALSE)
}