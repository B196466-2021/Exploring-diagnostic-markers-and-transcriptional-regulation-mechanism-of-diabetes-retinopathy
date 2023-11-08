#! /usr/bin/Rscript

usage<-function(spec){
	cat("This script is used for classified statistic and differential analysis.\n",
getopt(spec,usage=TRUE),
"Parameters:
	-f, --file	input file.
	-c, --class	column name of classification.
	-s, --score	column name of score.
 
 Options:
	-t, --type	plot type: boxplot, violin, beeswarm, point, violin,boxplot (default boxplot) or barplot (only for classified data).
	-m, --meth	differential test method: wilcox.test, t.test, between two group (default no test) or chisq.test, fisher.test (only for classified data).
	-l, --log2	log2 transform or not (default FALSE).
	-n, --num	diplay number of each group or not (default FALSE).
	-x, --xlab	label of X axis (defaut NULL).
	-y, --ylab	label of Y axis (default column name of score).
	-p, --pal	palette type of ggsci: npg, aaas, nejm, lancet, jama, jco, ucscgb, d3, locuszoom, igv, uchicago, startrek, tron, futurama, rickandmorty or simpsons (default npg).
	-r, --trim	trim the tails of the violins to the range of the data or not (default TRUE).
	-w, --bwi	box width of boxplot 0 ~ 1 (default 0.8).
	-a, --aste	translate differential p values into asterisk or not (default FALSE).
	-k, --angle	angle of text in X axis (defaut 0).
	-i, --title	plot title (default NULL).
	-g, --grid	remove the grid or not (default FALSE).
	-b, --lsize	the label and tite size of axis (defult 20).
	-e, --tsize	the text size of axis (defult 15).
	-d, --wid	the width of out plot (defult 7).
	-j, --hei	the height of out plot (defult 7).
	-o, --out	prefix of out files (default ./).
	-v, --version	display version information.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1)
}

##########################
#get options
library(getopt)
spec = matrix(c(
	'file','f',1,'character',
	'class','c',1,'character',
	'score','s',1,'character',
	'type','t',1,'character',
	'meth','m',1,'character',
	'log2','l',2,'logical',
	'num','n',2,'logical',
	'xlab','x',1,'character',
	'ylab','y',1,'character',
	'pal','p',1,'character',
	'trim','r',2,'logical',
	'bwi','w',1,'double',
	'aste','a',2,'logical',
	'title','i',1,'character',
	'grid','g',2,'logical',
	'lsize','b',1,'double',
	'tsize','e',1,'double',
	'wid','d',1,'double',
	'hei','j',1,'double',
	'angle','k',1,'double',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || is.null(opt$class) || is.null(opt$score) || !is.null(opt$help)) {usage(spec)}
if (is.null(opt$type)) {opt$type='boxplot'}
if (is.null(opt$log2)) {opt$log2=F}
if (is.null(opt$num)) {opt$num=F}
if (is.null(opt$ylab)) {opt$ylab=opt$score}
if (is.null(opt$pal)) {opt$pal='npg'}
if (is.null(opt$trim)) {opt$trim=T}
if (is.null(opt$bwi)) {opt$bwi=0.8}
if (is.null(opt$aste)) {opt$aste=F}
if (is.null(opt$grid)) {opt$grid=F}
if (is.null(opt$lsize)) {opt$lsize=20}
if (is.null(opt$tsize)) {opt$tsize=15}
if (is.null(opt$angle)) {opt$angle=0}
if (is.null(opt$wid)) {opt$wid=7}
if (is.null(opt$hei)) {opt$hei=7}
if (is.null(opt$out)) {opt$out='./'}
dir <- gsub("[^/]*$", '', opt$out, perl=T)
if(! dir.exists(dir) && dir != ''){
	dir.create(dir, recursive=T)
}
##########################
#define functions

##########################
#main programe
library(ggplot2)
library(beeswarm)
library(dplyr)
library("ggsci")
types <- strsplit(opt$type, ',')[[1]]
da <- read.table(opt$file, header=T, sep='\t', check.names=F, stringsAsFactors=F)
da2 <- da[!is.na(da[,opt$score]) & !is.na(da[,opt$class]), c(opt$score, opt$class)]
names(da2) <- c('score', 'class')
if(opt$log2){
	if(all((da2$score + 1) > 0)){
		da2$score <- log2(da2$score + 1)
	}else{
		stop('can not do log2 transform.')
	}
}
if(!'barplot' %in% types){
	pd <- beeswarm(score~class, data=da2, spacing=0.8)
	if(class(da2$class) %in% c("integer", "numeric")){
		pd$x.orig <- factor(pd$x.orig, levels=as.character(sort(unique(as.numeric(pd$x.orig)))))
	}
	p <- ggplot(data=pd)
	if('beeswarm' %in% types){
		p <- p + geom_point(aes(x=x, y=y, color=x.orig))
	}
	if('point' %in% types){
		p <- p + geom_jitter(aes(x=as.numeric(factor(x.orig)), y=y.orig, color=x.orig), width=0.3)
	}
	if('violin' %in% types){
		if(any(c('beeswarm', 'point') %in% types)){
			p <- p + geom_violin(aes(x=as.numeric(factor(x.orig)), y=y.orig, group=x.orig), fill=NA, trim=opt$trim)
		}else{
			p <- p + geom_violin(aes(x=as.numeric(factor(x.orig)), y=y.orig, fill=x.orig), trim=opt$trim)
			p <- p + get(paste('scale_fill_', opt$pal, sep=''))()
		}
	}
	if('boxplot' %in% types){
		if(any(c('beeswarm', 'point') %in% types)){
			p <- p + geom_boxplot(aes(x=as.numeric(factor(x.orig)), y=y.orig, group=x.orig), fill=NA, width=opt$bwi, outlier.shape=NA)
		}else if('violin' %in% types){
			p <- p + geom_boxplot(aes(x=as.numeric(factor(x.orig)), y=y.orig, group=x.orig), fill='white', width=opt$bwi, outlier.shape=NA)
		}else{
			p <- p + geom_boxplot(aes(x=as.numeric(factor(x.orig)), y=y.orig, fill=x.orig), linetype="dashed", width=opt$bwi, outlier.shape=NA)
			p <- p + stat_boxplot(aes(x=as.numeric(factor(x.orig)), y=y.orig, ymin=..lower..,ymax=..upper.., fill=x.orig) , width=opt$bwi, outlier.shape=NA)
			p <- p + stat_boxplot(geom = "errorbar", aes(x=as.numeric(factor(x.orig)), y=y.orig, ymin=..ymax.., group=x.orig), width=opt$bwi*0.6)
			p <- p + stat_boxplot(geom = "errorbar", aes(x=as.numeric(factor(x.orig)), y=y.orig, ymax=..ymin.., group=x.orig), width=opt$bwi*0.6)
			p <- p + get(paste('scale_fill_', opt$pal, sep=''))()
		}
	}
}else if('barplot' %in% types){
	pd <- da2 %>% group_by(class, score) %>% summarise(count=length(score)) %>% group_by(class) %>% mutate(total=sum(count), ratio=count/total) %>% data.frame()
	pd$x.orig <- pd$class
	pd$x <- as.numeric(factor(pd$x.orig))
	p <- ggplot(pd, aes(x=x, y=ratio, fill=score)) + geom_col()
}

gp.stat <- ggplot_build(p)$data
pd.md <-  data.frame()
for(gs in gp.stat){
	sub.d <- gs[, names(gs) %in% c('x', 'y', 'ymin', 'ymax')]
    md <- data.frame(x=rep(round(sub.d$x), ncol(sub.d)-1), y=as.numeric(as.matrix(sub.d[grep('^y', names(sub.d))])))
	pd.md <- rbind(pd.md, md)
	pd.md <- pd.md[!is.na(pd.md$y),]
}
pd.stat <- pd.md %>% group_by(x) %>% summarise(ymin=min(y), ymax=max(y)) %>% data.frame()
labe.caption <- NULL
if(!is.null(opt$meth)){
	cbn <- utils::combn(sort(unique(as.numeric(factor(pd$x.orig)))), m=2)
	p.dat <- data.frame()
	for(i in 1:ncol(cbn)){
		da.sub <- pd[as.numeric(factor(pd$x.orig)) %in% cbn[, i], ]
		if(!'barplot' %in% types){
			dre <- get(opt$meth)(formula=y.orig ~ x.orig, data=da.sub)
		}else if('barplot' %in% types){
			dre <- get(opt$meth)(table(da2))
		}
		pos <- cbn[, i]
		p.sub <- data.frame(x1=min(pos), x=mean(pos), x2=max(pos), pv=dre$p.value, stringsAsFactors=F)
		p.dat <- rbind(p.dat, p.sub)
	}
	da.label <- data.frame(t(apply(p.dat, 1, function(x){x=c(x, max(pd.stat[x[1]:x[3], 'ymax']) + (max(pd.stat$ymax)-min(pd.stat$ymin))*0.05)})))
	names(da.label)[ncol(da.label)] <- 'y'
	da.label$label <- format(da.label$pv, digits=2, scientific=T)
	if(opt$aste){
		da.label <- da.label[da.label$pv < 0.05,]
		da.label$label[da.label$pv < 0.05 & da.label$pv >= 0.01] <- '*'
		da.label$label[da.label$pv < 0.01 & da.label$pv >= 0.001] <- '**'
		da.label$label[da.label$pv < 0.001] <- '***'
		labe.size <- 8
		labe.caption <- expression(paste('*: ', 0.01 <= p, '< 0.05;  **: ', 0.001 <= p, '< 0.01;  ***: p < 0.001'))
	}else{
		labe.size <- 5
	}
	if(nrow(da.label) >0){
		da.label$order <- as.numeric(factor(da.label$x2 - da.label$x1))
		da.label <- da.label[order(da.label$order, da.label$x1, da.label$x2), ]
		da.label$coef <- 0
		if(nrow(da.label) > 1){
			for(i in 2:nrow(da.label)){
				cors <- (da.label$x1[1:i-1] < da.label$x1[i] & da.label$x2[1:i-1] > da.label$x1[i]) | (da.label$x1[1:i-1] < da.label$x2[i] & da.label$x2[1:i-1] > da.label$x2[i]) | 
				(da.label$x1[1:i-1] < da.label$x2[i] & da.label$x1[1:i-1] > da.label$x1[i]) | (da.label$x2[1:i-1] < da.label$x2[i] & da.label$x2[1:i-1] > da.label$x1[i])
				if(any(cors)){
					da.label$coef[i] <- max(da.label$coef[1:i-1][cors]) + 1
				}
			}
		}
		da.label$y.label <- da.label$y + (max(pd.stat$ymax)-min(pd.stat$ymin))*0.08*da.label$coef
		da.label$y2 <- da.label$y.label - (max(pd.stat$ymax)-min(pd.stat$ymin))*0.03
		da.label$y1 <- da.label$y.label - (max(pd.stat$ymax)-min(pd.stat$ymin))*0.04
		sld <- data.frame(x1=c(da.label$x1+0.05, da.label$x2-0.05), y1= rep(da.label$y1, 2), x2=c(da.label$x1+0.05, da.label$x2-0.05), y2=rep(da.label$y2, 2))
		p <- p + geom_text(data=da.label, aes(x=x, y=y.label, label=label), size=labe.size, inherit.aes=F)
		p <- p + geom_segment(aes(x = x1+0.05, y = y2, xend = x2-0.05, yend = y2), data = da.label, inherit.aes=F)
		p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data=sld, inherit.aes=F)
		ylimt <- c(min(pd.stat$ymin), max(c(da.label$y.label + (max(pd.stat$ymax)-min(pd.stat$ymin))*0.01, pd.stat$ymax)))
	}else{
		ylimt <- c(min(pd.stat$ymin), max(pd.stat$ymax))
		labe.caption <- NULL
	}
}else{
	ylimt <- c(min(pd.stat$ymin), max(pd.stat$ymax))
}
xls <- pd[,c('x', 'x.orig')]
xls$x <- as.numeric(factor(xls$x.orig))
xls$lab <- xls$x.orig
if(opt$num){
	xls <- xls %>% group_by(x) %>% mutate(num=length(x.orig), lab=paste(x.orig, '\n(', num, ')', sep='')) %>% data.frame(stringsAsFactors=F)
}
xlsb <- unique(xls)
p <- p + scale_y_continuous(name=opt$ylab, expand=c(0.03, 0)) + coord_cartesian(ylim=ylimt)
p <- p + scale_x_continuous(name=opt$xlab, breaks=xlsb$x, labels=xlsb$lab)
p <- p + get(paste('scale_fill_', opt$pal, sep=''))()
p <- p + labs(title=opt$title, caption=labe.caption) + theme_bw()
if(opt$grid){
	p <- p + theme(panel.grid=element_blank())
}
if(!'barplot' %in% types){
	p <- p + theme(axis.title=element_text(size=opt$lsize), axis.text=element_text(size=opt$tsize, color='black'), 
				   plot.title=element_text(hjust=0.5, size=opt$lsize), legend.position='none')
}else if('barplot' %in% types){
	p <- p + theme(axis.title=element_text(size=opt$lsize), axis.text=element_text(size=opt$tsize, color='black'),
				   plot.title=element_text(hjust=0.5, size=opt$lsize), legend.title=element_blank(), legend.text=element_text(size=opt$tsize))
}
if(opt$angle != 0){
	p <- p + theme(axis.text.x=element_text(size=opt$tsize, color='black', angle=opt$angle, hjust=1))
}
ggsave(paste(opt$out, gsub(',', '_', opt$type), '.png', sep=''), plot=p, width=opt$wid, height=opt$hei, units='in', dpi=600)
ggsave(paste(opt$out, gsub(',', '_', opt$type), '.pdf', sep=''), plot=p, width=opt$wid, height=opt$hei, units='in')
if(file.exists('Rplots.pdf')){
	file.remove('Rplots.pdf')
}
