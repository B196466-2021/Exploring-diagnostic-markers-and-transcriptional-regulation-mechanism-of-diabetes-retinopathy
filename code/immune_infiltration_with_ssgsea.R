#! /usr/bin/Rscript

usage<-function(spec){
	cat("This script is used to estimate immune infiltration with ssGSEA method\n",
getopt(spec,usage=TRUE),
"Options:
	-e, --exp	expression file with non-negative data (counts or log2 transformed normalized expression).
	-d, --db	database dataset (immune cell<tab>NA<table>gene1<table>gene2...).
	-p, --pheno	group information file of samples.
	-m, --method	gsva method gsva/ssgsea/zscore/plage (default ssgsea).
	-l, --log2	convert expression value by log2(x+1) or not (default FALSE).
	-t, --thread	threads number (default 5)
	-o, --out	output directory (deafault ./).
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
	'db','d',1,'character',
	'pheno','p',1,'character',
	'thread','t',1,'integer',
	'method','m',1,'character',
	'log2','l',2,'logical',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$exp) || is.null(opt$db) || !is.null(opt$help)) {usage(spec)}

if (is.null(opt$out)) {opt$out='./'}
if (is.null(opt$thread)) {opt$thread=5}
if (is.null(opt$method)) {opt$method='ssgsea'}
if (is.null(opt$log2)) {opt$log2=F}
dir <- gsub("[^/]*$", '', opt$out, perl=T)
if(! dir.exists(dir) && dir != ''){
	dir.create(dir, recursive=T)
}
##########################
#define functions

my_fun<-function(){

}

##########################
#main programe
library(genefilter)
library(GSVA)
library(pheatmap)
library(limma)
library(ggplot2)
library(edgeR)

db <- readLines(opt$db)
dbl <- list()
for(i in db){
	ic <- strsplit(i, split='\t')[[1]]
	dbl[[ic[1]]] <- ic[-c(1, 2)]
}
exp <- read.table(opt$exp, header=T, sep='\t', check.names=F, stringsAsFactors=F)
if(!is.null(opt$pheno) && file.exists(opt$pheno)){
	ga <- read.table(opt$pheno, header=T, sep='\t', check.names=F, stringsAsFactors=F)
	group_name <- names(ga)[2]
	names(ga)[2] <- 'Group'
	exp <- exp[, names(exp) %in% c(names(exp)[1], ga$Sample)]
}

if(any(duplicated(exp[,1]))){
	exp <- aggregate(exp[,-1], by=list(exp[,1]), mean)
	rownames(exp) <- exp[,1]
	exp <- exp[,-1]
}else{
	rownames(exp) <- exp[,1]
	exp <- exp[,-1]
}
expt <- as.matrix(exp)
if(opt$log2){
	expt <- log2(expt + 1)
}
expt <- expt[! apply(expt, 1, function(x){all(x==0)}),]
if(length(which(expt%%1==0))/length(expt) > 0.8){
	expn <- log2(cpm(expt) + 1)
	expn <- expn[genefilter(expn, filterfun(kOverA(round(ncol(expn)*0.05), 0))), ]
}else if(max(expt[!is.na(expt)]) > 25){
	expn <- log2(expt + 1)
	expn <- expn[genefilter(expn, filterfun(kOverA(round(ncol(expn)*0.05), 0))), ]
}else{
	expn <- expt
	expn <- expn[genefilter(expn, filterfun(kOverA(round(ncol(expn)*0.05), 0))), ]
}

if(length(which(expn%%1==0))/(nrow(expn)*ncol(expn)) > 0.8){
	kcdf <- 'Poisson'
}else{
	kcdf <- 'Gaussian'
}
gsva.re <- gsva(expn, dbl, method=opt$method, kcdf=kcdf, abs.ranking=TRUE, parallel.sz=opt$thread, verbose=FALSE)
gsva.or <- gsva.re[order(apply(gsva.re, 1, mean), decreasing=T), ]
if(class(gsva.or)=="numeric"){
	od <- data.frame(names(gsva.or), gsva.or, check.names=F)
	names(od) <- c('id', names(dbl)[1])
}else{
	od <- data.frame(id=rownames(gsva.or), gsva.or, check.names=F)
}
write.table(od, file=paste(opt$out, '/' , opt$method, '_result.txt', sep=''), sep='\t', quote=F, row.names=F)
##correlation
if(nrow(gsva.re) >1){
	library(psych)
	cor.re <- corr.test(t(gsva.or), adjust='none')
	cor.pd <- data.frame(x=rep(rownames(cor.re$r), ncol(cor.re$r)), y=rep(colnames(cor.re$r), each=nrow(cor.re$r)), 
						 r=as.numeric(cor.re$r), p=round(as.numeric(cor.re$p), 2), stringsAsFactors=F)
	cor.pd$x <- factor(cor.pd$x, levels=rownames(cor.re$r))
	cor.pd$y <- factor(cor.pd$y, levels=rev(colnames(cor.re$r)))
	cor.pd$p <- as.character(cor.pd$p)
	cor.pd$p[is.na(cor.pd$p)] <- 'NA'

	p <- ggplot(cor.pd, aes(x, y, fill=r)) + geom_tile(colour='grey50')
	p <- p + scale_fill_gradientn(colors=c('blue','white','red'), na.value='gray70', n.breaks=3, 
								  guide=guide_colourbar(title='R', title.position='top', label.position='bottom', direction='horizontal', title.hjust=0.5, ticks=F),
								  values = scales::rescale(c(min(cor.pd$r[!is.na(cor.pd$r)]), 0, max(cor.pd$r[!is.na(cor.pd$r)]))))
	p <- p + scale_y_discrete(position='right', expand=c(0, 0)) + scale_x_discrete(, expand=c(0, 0))
	p <- p + geom_text(aes(label=p), size=3) 
	p <- p + theme(axis.title=element_blank(), axis.ticks.length=unit(0, 'mm'), axis.text=element_text(size=12, colour='black'), 
				   axis.text.x=element_text(angle=90, hjust=1), legend.position=c(1.02, -0.02), legend.justification=c(0, 1), legend.background=element_blank())
	nw <- ifelse(nrow(gsva.re)>15, nrow(gsva.re), 15)
	#ggsave(paste(opt$out, '/' , opt$method, '_correlation_heatmap.pdf', sep=''), plot=p,  width=nw*0.4, height=nw*0.4, units='in')
	#ggsave(paste(opt$out, '/' , opt$method, '_correlation_heatmap.png', sep=''), plot=p,  width=nw*0.4, height=nw*0.4, units='in', dpi=600)
}
##diff box plot
if(nrow(gsva.re) > 15){
	show_coln <- F
}else{
	show_coln <- T
}
if(exists('ga') && length(unique(ga$Group))>1){
	gas <- ga[ga$Sample %in% colnames(gsva.or), ]
	gaso <- gas[order(gas$Group),]
	gsva.sub <- gsva.or[, as.character(gaso$Sample)]
	anno <- data.frame(Group=gaso$Group)
	rownames(anno) <- gaso$Sample
	legd <- T
	#diff box plot
	if(length(unique(gaso$Group))==2){
		pvalue <- apply(gsva.sub, 1, function(x){sub.da=data.frame(sc=x, group=gaso$Group); wtr=wilcox.test(sc ~ group, data=sub.da); return(wtr$p.value)})
	}else{
		pvalue <- apply(gsva.sub, 1, function(x){sub.da=data.frame(sc=x, group=gaso$Group); wtr=kruskal.test(sc ~ group, data=sub.da); return(wtr$p.value)})
	}
	ptxt <- data.frame(x1=names(pvalue), y1=max(gsva.sub[!is.na(gsva.sub)]), pval=pvalue)
	ptxt$pval[pvalue < 0.001] <- '***'
	ptxt$pval[pvalue >= 0.001 & pvalue < 0.01] <- '**'
	ptxt$pval[pvalue >= 0.01 & pvalue < 0.05] <- '*'
	ptxt$pval[pvalue >= 0.05] <- ''
	pre <- data.frame(x=rep(rownames(gsva.sub), ncol(gsva.sub)), y=as.vector(gsva.sub), Group=rep(gaso$Group, each=nrow(gsva.sub)), stringsAsFactors=F)
}else{
	anno <- NA
	gsva.sub <- gsva.or
	legd <- F
	pre <- data.frame(x=rep(rownames(gsva.sub), ncol(gsva.sub)), y=as.vector(gsva.sub), Group=1, stringsAsFactors=F)
}

p <- ggplot(pre, aes(x, y, fill=Group)) + geom_boxplot(show.legend=legd, outlier.shape=NA)
if(exists('pvalue')){
	p <- p + geom_text(aes(x=x1, y=y1, label=pval), data=ptxt, inherit.aes=F)
	p <- p + guides(fill=guide_legend(title=group_name))
}
p <- p + theme_bw() +  labs(x=NULL, y='Enrichment Scores')
p <- p + theme(axis.title=element_text(size=20), axis.text=element_text(size=15, colour='black'), axis.text.x=element_text(angle=45, hjust=1), 
			   legend.text=element_text(size=15), legend.title=element_blank())
w <- ifelse(nrow(gsva.re) > 3, nrow(gsva.re) + 2, 5)
ggsave(paste(opt$out, '/', opt$method, '_boxplot.pdf', sep=''), plot=p, width=w, height=5)
ggsave(paste(opt$out, '/', opt$method, '_boxplot.png', sep=''), plot=p, width=w, height=5, dpi=600)

##pheatmap
if(nrow(gsva.re) >1){
	wd <- 7 + max(nchar(rownames(gsva.sub)))*0.1
	p <- pheatmap(gsva.sub, scale='none', cluster_rows=F, cluster_cols=F, show_colnames=show_coln, annotation_col=anno, silent=T)
	ggsave(paste(opt$out, '/', opt$method, "_heatmap.png",sep=""), plot=p, width=wd, height=7, units='in', dpi=600)
	ggsave(paste(opt$out, '/', opt$method, "_heatmap.pdf",sep=""), plot=p, width=wd, height=7, units='in')
}

#oave.image()
