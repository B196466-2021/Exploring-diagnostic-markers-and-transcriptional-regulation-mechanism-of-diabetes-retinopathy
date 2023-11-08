#! /usr/bin/Rscript

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
