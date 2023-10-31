#! /usr/bin/Rscript

usage<-function(){
	cat("This script is used for co_expression analysis with WGCNA.\n",
getopt(spec,usage=TRUE),
"Options:
	-e, --exp	expression file.
			ID	sample1	sample2	sample3	...
			gene1	31	37	31	...
			gene2	23	34	65	...
	-p, --phen	phenotype or clinical traits file (1 is the trait, 0 is not the trait).
			Sanple	status
			sample1	1
			sample2	1
			sample3	0
	-r, --corr	correlation threshold, must be -1 ~ 1 (default 0.8).
	-c, --count	gene count selected (default 5000 genes of SD).
	-n, --norm	normalize or not by limma (default FALSE).
	-t, --thr	thread will be used (at least 2, default 2).
	-o, --out	out directory (default ./).
	-v, --version	display version information.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1)
}

##########################
#get options

library(getopt)
spec = matrix(c(
'exp','e',1,'character',
'phen','p',1,'character',
'corr','r',1,'double',
'thr','t',1,'integer',
'count','c',1,'integer',
'norm','n',0,'logical',
'out','o',1,'character',
'version','v',0,'logical',
'help','h',0,'logical'
),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$exp) || is.null(opt$phen) || !is.null(opt$help)) {usage()}

if (is.null(opt$corr)) {opt$corr=0.8}
if (is.null(opt$thr)) {opt$thr=2}
if (is.null(opt$out)) {opt$out='./'}
if (is.null(opt$count)) {opt$count=5000}
if (is.null(opt$norm)) {opt$norm=FALSE}
dir <- gsub("[^/]*$", '', opt$out, perl=T)
if(! dir.exists(dir) && dir != ''){
	dir.create(dir, recursive=T)
}

##########################
#main programe
library(limma)
library(WGCNA)
library(ggplot2)
library(cowplot)
library(psych)
options(stringsAsFactors = FALSE)
enableWGCNAThreads(opt$thr)
source('/home/dev/DEV-wangdh/R/labeledHeatmap.R')
source('/home/dev/DEV-wangdh/R/heatmapWithLegend.R')

## Loading datas
exp<-read.table(opt$exp, header=T, sep='\t', check.names=F, stringsAsFactors=F)
if(any(duplicated(exp[, 1]))){
	exp <- aggregate(exp[, -1], by=list(exp[, 1]), function(x)(mean(x[!is.na(x)])))
}
rownames(exp) <- exp[,1]
exp <- exp[,-1]
phen <- as.matrix(read.table(opt$phen, header=T, sep='\t', check.names=F, stringsAsFactors=F, row.names=1))
exps<-exp[, names(exp) %in% rownames(phen)]
phe <- data.frame(phen[names(exps),], stringsAsFactors=F)
names(phe) <- colnames(phen)
rownames(phe) <- names(exps)
if(length(which(as.matrix(exps)%%1==0))/length(as.matrix(exps)) > 0.8){
	exps<-voom(exps)$E
}
if(opt$norm){
	print('Normalize expression datas ...')
	expn<-normalizeBetweenArrays(exps, method='quantile')
}else{
	expn <- exps
}
expt<-as.data.frame(t(expn))
if(ncol(expt) > opt$count){
	gsdr <- apply(expt, 2, function(x){sd(x, na.rm=T)})
	expt <- expt[, names(head(sort(gsdr, decreasing=T), opt$count))]
}
## Checking data for excessive missing values
gsg = goodSamplesGenes(expt, verbose = 3)
if(!gsg$allOK){
	if(sum(!gsg$goodGenes)>0){
		printFlush(paste("Removing genes:", paste(names(expt)[!gsg$goodGenes], collapse = ", ")))
	}
	if(sum(!gsg$goodSamples)>0){
		printFlush(paste("Removing samples:", paste(rownames(expt)[!gsg$goodSamples], collapse = ", ")))
	}
	expt<-expt[gsg$goodSamples, gsg$goodGenes]
}
## Processing clinical trait data
phe2<-data.frame(phe[rownames(phe) %in% rownames(expt), ], check.names=F)
names(phe2)<-names(phe)
head(phe2)
rownames(phe2) <- rownames(phe)[rownames(phe) %in% rownames(expt)]
phe2 <- data.frame(phe2[rownames(expt), ], check.names=F)
names(phe2) <- names(phe)
rownames(phe2) <- rownames(expt)
## Sample clust tree
sampletree = hclust(dist(expt), method = "average")
str(phe2)
traitcolors = numbers2colors(phe2, signed = FALSE)
tsz <- 50/nrow(phe2)
if(tsz > 1){
	tsz <- 1
}
pdf(paste(opt$out, '/Dendrogram_sample_pheno_traits.pdf', sep=''), width = 10)
plotDendroAndColors(sampletree, traitcolors, groupLabels = names(phe2), main = "Sample dendrogram and trait heatmap", cex.dendroLabels=tsz, cex.colorLabels=0.7)
dev.off()
png(paste(opt$out, '/Dendrogram_sample_pheno_traits.png', sep=''), width = 10, height = 7, units = "in", res=600)
plotDendroAndColors(sampletree, traitcolors, groupLabels = names(phe2), main = "Sample dendrogram and trait heatmap", cex.dendroLabels=tsz, cex.colorLabels=0.7)
dev.off()

## Automatic construction of the gene network and identication of modules
powers<-c(c(1:10), seq(from = 12, to=30, by=2))
sft<-pickSoftThreshold(expt, powerVector = powers, verbose = 5)

# determining the power
power_candidate_list=-sign(sft$fitIndices[,3])*sft$fitIndices[,2]
print('Correlation value of each powers:')
print(sort(power_candidate_list, decreasing=T))

if(any(power_candidate_list>opt$corr)){
	power_select<-min(powers[power_candidate_list>opt$corr])
	corr_select<-min(power_candidate_list[power_candidate_list>opt$corr])
}else{
	power_select<-powers[order(power_candidate_list, decreasing=T)[1]]
	corr_select<-max(power_candidate_list)
}

print(paste('Selected power is: ', power_select, ' correlation value is: ', corr_select, ''))
psd1<-data.frame(x=sft$fitIndices[,1], y=-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
psd2<-data.frame(x=sft$fitIndices[,1], y=sft$fitIndices[,5])
p <- ggplot(psd1, aes(x, y)) + geom_point()
p <- p + geom_hline(yintercept = corr_select, colour = "red")
p <- p + labs(x='Soft Threshold (power)', y=expression('Scale Free Topology Model Fit,signed R'^{2}), title='Scale independence')
p <- p + theme(plot.title=element_text(hjust = 0.5), panel.background = element_rect(fill = "white", colour = "black"))
p <- p + geom_vline(xintercept = power_select, linetype=2, colour = "red")
p1 <- p + annotate('text', x=c(power_select, 3), y=c(0.5, corr_select + 0.01), label=c(paste('Power = ', power_select, sep=''), as.expression(substitute(italic('R')^2~"="~r2, list(r2=round(corr_select, 3))))), colour = "red")

p <- ggplot(psd2, aes(x, y)) + geom_point()
p <- p + labs(x='Soft Threshold (power)', y='Mean Connectivity', title='Mean connectivity')
p2 <- p + theme(plot.title=element_text(hjust = 0.5), panel.background = element_rect(fill = "white", colour = "black"))

pdf(file=paste(opt$out, '/Power_select.pdf', sep=''), width=9)
plot_grid(p1, p2, ncol=2)
dev.off()

png(file=paste(opt$out, '/Power_select.png', sep=''), height=7, width=9, res=600, units='in')
plot_grid(p1, p2, ncol=2)
dev.off()

## Step-by-step network construction
adjacency<-adjacency(expt, power=power_select)

TOM<-TOMsimilarity(adjacency)

dissTOM<-1-TOM

geneTree<-hclust(as.dist(dissTOM), method="average")

dynamicMods<-cutreeDynamic(dendro=geneTree, distM=dissTOM,
							deepSplit=2, pamRespectsDendro=FALSE,
							minClusterSize=30)

dynamicColors<-labels2colors(dynamicMods)

#Merging modules
MEList<-moduleEigengenes(expt, colors=dynamicColors)
MEs<-MEList$eigengenes
MEDiss<-1-cor(MEs)
METree<-hclust(as.dist(MEDiss), method="average")
cutv<-0.2
mergemodule<-mergeCloseModules(expt, dynamicColors, verbose = 3, cutHeight = cutv)
while(ncol(mergemodule$newMEs) > 20){
	cutv <- cutv + 0.1
	mergemodule<-mergeCloseModules(expt, dynamicColors, verbose = 3, cutHeight = cutv)
}
mergedColors<-mergemodule$colors
mergedMEs<-mergemodule$newMEs
nMEs<-mergemodule$newMEs

mogn<-as.data.frame(table(mergedColors))
names(mogn)<-c('Module', 'Gene_counts')
write.table(mogn, paste(opt$out, '/Gene_counts_for_each_module.txt', sep=''), quote=F, sep='\t', row.names=F)

# Plot the dendrogram and the module colors underneath
pdf(paste(opt$out, '/Dendrogram_coexpression_gene_clustering_result.pdf', sep=''), width=14)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
					c("Module", "Merged module"),
					dendroLabels = FALSE, hang = 0.03,
					addGuide = TRUE, guideHang = 0.05,
					main="Gene dendrogram and module colors")

dev.off()
png(paste(opt$out, '/Dendrogram_coexpression_gene_clustering_result.png', sep=''), width = 14, height = 7, units = "in", res=600)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
					c("Module", "Merged module"),
					dendroLabels = FALSE, hang = 0.03,
					addGuide=TRUE, guideHang=0.05,
					main="Gene dendrogram and module colors")

dev.off()

## Visualizing the gene network
plotTOM<-dissTOM^7
diag(plotTOM)<-NA
pdf(paste(opt$out, '/Heatmap_gene_gene_correlation.pdf', sep=''))
TOMplot(plotTOM, geneTree, mergedColors, main = "Network heatmap plot, all genes")
dev.off()
png(paste(opt$out, '/Heatmap_gene_gene_correlation.png', sep=''), width = 7, height = 7, units = "in", res=600)
TOMplot(plotTOM, geneTree, mergedColors, main = "Network heatmap plot, all genes")
dev.off()

## Visualizing the network of eigengenes
MET<-orderMEs(cbind(nMEs, phe2))
pdf(paste(opt$out, '/Relationships_among_the_eigengenes_and_the_trait.pdf', sep=''))
plotEigengeneNetworks(MET, "", marDendro = c(0,5.5,1,3), marHeatmap = c(5,6,1,1), cex.lab = 0.8, xLabelsAngle= 30, excludeGrey=F)
dev.off()
png(paste(opt$out, '/Relationships_among_the_eigengenes_and_the_trait.png', sep=''), width = 7, height = 7, units = "in", res=600)
plotEigengeneNetworks(MET, "", marDendro = c(0,5.5,1,3), marHeatmap = c(5,6,1,1), cex.lab = 0.8, xLabelsAngle= 30, excludeGrey=F)
dev.off()

## Exporting network data
probes<-names(expt)
for(i in unique(mergedColors)){
	inModule<-is.finite(match(mergedColors, i))
	modProbes<-probes[inModule]
	modTOM<-TOM[inModule, inModule]
	cdir<-paste(opt$out, '/cytoscape_network/', i, sep='')
	dir.create(cdir, recursive=TRUE)
	cyt<-exportNetworkToCytoscape(modTOM, edgeFile=paste(cdir, "/", i, "_edges.txt", sep=""),
								  nodeFile=paste(cdir, "/", i, "_nodes.txt", sep=""),
								  weighted=TRUE, threshold=0, nodeNames=modProbes,
								  nodeAttr=mergedColors[inModule])
}
gme<-data.frame(ID=probes, Module=mergedColors, t(expt), check.names=F)
gme<-gme[order(gme$Module),]
write.table(gme, paste(opt$out, '/Gene_clustering_result.txt', sep=''), quote=F, sep='\t', row.names=F)
## Relating modules to external clinical traits
file.copy(from='/home/share/Readme/Readme_WGCNA.pdf', to=paste(opt$out, '/', sep=''))

moduleTraitCor<-cor(orderMEs(nMEs), phe2, use = "p")
colnames(moduleTraitCor) <- paste(colnames(moduleTraitCor), '(cor)', sep='')
moduleTraitPvalue<-corPvalueStudent(moduleTraitCor, nrow(phe2))
colnames(moduleTraitPvalue) <- paste(colnames(moduleTraitCor), '(pval)', sep='')
copv<-merge(moduleTraitCor, moduleTraitPvalue, by=0)
names(copv)[1]<-'Module'
write.table(copv, paste(opt$out, '/Module_trait_correlation_significance.txt', sep=''), quote=F, sep='\t', row.names=F)

pla <- moduleTraitPvalue
pla[moduleTraitPvalue < 0.01] <- '\n***'
pla[moduleTraitPvalue >= 0.01 & moduleTraitPvalue < 0.05] <- '\n**'
pla[moduleTraitPvalue >= 0.05 & moduleTraitPvalue < 0.1] <- '\n*'
pla[moduleTraitPvalue >= 0.1] <- ''

textMatrix = paste(signif(moduleTraitCor, 2), pla, sep = "")

#dim(textMatrix)<- dim(moduleTraitCor)
if(ncol(moduleTraitCor) >10){
	w <- 7*ncol(moduleTraitCor)/15
}else{
	w <- 7
}
if(row(moduleTraitCor) >15){
	h <- 7*nrow(moduleTraitCor)/15
}else{
	h <- 7
}
#save.image(file='co_exp.RData')
##plot correlation between module and trait
pd1 <- data.frame(me=rownames(moduleTraitCor), color=gsub('^ME', '', rownames(moduleTraitCor)), stringsAsFactors=F) 
pd1$me <- factor(pd1$me, levels=rev(rownames(moduleTraitCor)))
pd2 <- data.frame(me=rep(rownames(moduleTraitCor), ncol(moduleTraitCor)), corr=as.numeric(moduleTraitCor), tarit=rep(gsub('\\(cor\\)', '', colnames(moduleTraitCor)), each=nrow(moduleTraitCor)), lab=textMatrix, stringsAsFactors=F)
pd2$me <- factor(pd2$me, levels=rev(rownames(moduleTraitCor)))
pd2$tarit <- factor(pd2$tarit, levels=unique(pd2$tarit)[order(nchar(unique(pd2$tarit)))])

p1 <- ggplot(pd1, aes(x='a', y=me, fill=color)) + geom_tile() + scale_fill_identity()
p1 <- p1 + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))
p1 <- p1 + theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.text.y=element_text(size=15, color='black'), axis.ticks=element_blank(), panel.border=element_blank(), plot.margin=margin(t=3, r=0, b=2, l=2, unit='mm'), plot.background=element_rect(fill='white'))

p2 <- ggplot(pd2, aes(x=tarit, y=me, fill=corr, label=lab)) + geom_tile() + geom_text(size=3)
p2 <- p2 + scale_fill_gradient2(low ="#0D8CFF", mid ="white", high ="#FF3300", limits=c(-1, 1), guide=guide_colourbar(title=NULL, barwidth = 1, barheight = 30))
p2 <- p2 + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))
p2 <- p2 + theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(size=15, color='black', angle=30, hjust=1), axis.ticks=element_blank(), plot.margin=margin(t=2, r=2, b=2, l=0.5, unit='mm'), legend.text=element_text(size=15), panel.border=element_rect(fill=NA, colour="black"), plot.background=element_rect(fill='white'))
p <- plot_grid(p1, p2, ncol=2, align='h', rel_widths=c(2, 5))
ggsave(paste(opt$out, '/Heatmap_and_barplot_module_trait_correlation.pdf', sep=''), plot=p, width=7, height=7, units='in')
ggsave(paste(opt$out, '/Heatmap_and_barplot_module_trait_correlation.png', sep=''), plot=p, width=7, height=7, units='in', dpi=600)

#nodule select
ADJ1 <- abs(cor(expt,use="p"))^sft$powerEstimate
Alldegrees1 <- intramodularConnectivity(ADJ1, mergedColors)
datKME <- signedKME(expt, nMEs, outputColumnName="MM.")

for(i in names(phe2)){
	GS <- as.numeric(cor(phe2[, i], expt, use="p"))
	GeneSignificance <- abs(GS)
	ModuleSignificance <- tapply(GeneSignificance, mergedColors, mean, na.rm=T)
	png(paste(opt$out, '/trait_', i, '_correlation_cignificance_barplot.png', sep=''), width = 9, height = 7, units = "in", res=600)
	par(mar=c(10, 5, 4, 2), cex.axis=1.5, cex.lab=1.5)
	plotModuleSignificance(GeneSignificance, mergedColors, las=2)
	dev.off()
	
	pdf(paste(opt$out, '/trait_', i, '_correlation_cignificance_barplot.pdf', sep=''), width = 9, height = 7)
	par(mar=c(10, 5, 4, 2), cex.axis=1.5, cex.lab=1.5)
	plotModuleSignificance(GeneSignificance, mergedColors, las=2)
	dev.off()
}

# correlation between genes and modules
gmd <- data.frame()
for(n in unique(gme$Module)){
	mexp <- t(gme[gme$Module==n, rownames(MEs)])
	mmes <- matrix(MEs[, paste('ME', n, sep='')], ncol=1)
	rownames(mmes) <- rownames(MEs)
	colnames(mmes) <- n
	corre <- corr.test(x=mexp, y=mmes)
	cored <- data.frame(Gene=rownames(corre$r), corr=corre$r[,1], p.value=corre$p[,1], p.adj=corre$p.adj[,1], Module=paste('ME', n, sep=''), stringsAsFactors=F, check.names=F)
	gmd <- rbind(gmd, cored)
}
write.table(gmd, paste(opt$out, '/gene_module_correlation.txt', sep=''), sep='\t', quote=F, row.names=F)
# correlation between genes and traits
core <- corr.test(x=expt, y=phe2)
ord <- data.frame(Gene=rownames(core$r), core$r, stringsAsFactors=F, check.names=F)
orp <- data.frame(Gene=rownames(core$p), core$p, stringsAsFactors=F, check.names=F)
write.table(ord, paste(opt$out, '/gene_trait_correlation_corr.txt', sep=''), sep='\t', quote=F, row.names=F)
write.table(orp, paste(opt$out, '/gene_trait_correlation_p_vlaue.txt', sep=''), sep='\t', quote=F, row.names=F)
q(save ='yes')
