#! /usr/bin/Rscript

##########################
#main programe
library(ggplot2)
library(dplyr)
library(ggsci)
colors <- strsplit(opt$color, ',')[[1]]
if(length(colors)==1){
	colors <- get(paste('pal_', colors, sep=''))()(7)
}
dd <- read.table(opt$data, header=T, sep='\t', check.names=F, row.names=1)
da <- as.matrix(dd)
gr <- read.table(opt$group, header=T, sep='\t', check.names=F, row.names=1, stringsAsFactors=F)
gr$Group <- factor(gr$Group, levels=unique(gr$Group))
if(length(which(colnames(da) %in% rownames(gr)))/nrow(gr) < 0.01){
	da <- t(da)
}

dtr <- data.frame()
dtd <- merge(t(da), gr, by=0)
da <- da[, dtd[,1]]
if(!any(c('p_value', 'p.value', 'p') %in% colnames(dd))){
	method <- get(opt$method)
	for(n in rownames(da)){
		dt <- method(formula=as.formula(paste('`', n, '` ~ Group', sep='')), data=dtd)
		std <- data.frame(drug=n, p.value=dt$p.value)
		dtr <- rbind(dtr, std)
	}
}else{
	dtr <- data.frame(drug=rownames(dd), p.value=dd[, which(colnames(dd) %in% c('p_value', 'p.value', 'p'))[1]])
}
md <- aggregate(dtd[, rownames(da)], by=list(dtd$Group), function(x)(mean(x[!is.na(x)])))
mdt <- t(md[,-1])
colnames(mdt) <- md[,1]
od <- merge(dtr, mdt, by.x=1, by.y=0)
names(od)[1] <- 'ID'
write.table(od, paste(opt$out, 'diff_test_result.txt', sep=''), sep='\t', quote=F, row.names=F)

if(opt$jdt){quit(save = "no")}

dtr$lab <- format(dtr$p.value, digits=2, scientific=T)
dtr$x <- as.numeric(factor(dtr$drug))
if(opt$switch){
	dtr$lab <- 'ns'
	dtr$lab[dtr$p.value < 0.001] <- '***'
	dtr$lab[dtr$p.value < 0.01 & dtr$p.value >= 0.001] <- '**'
	dtr$lab[dtr$p.value < 0.05 & dtr$p.value >= 0.01] <- '*'
	labe.caption <- expression(paste('*: ', 0.01 <= p, '< 0.05;  **: ', 0.001 <= p, '< 0.01;  ***: p < 0.001'))
	#dtr <- dtr[dtr$p.value < 0.1,]
}
dtr <- dtr[order(dtr$x), ]
pd <- data.frame(drug=rep(rownames(da), ncol(da)), sample=rep(colnames(da), each=nrow(da)), score=as.numeric(da), stringsAsFactors=F)
pd$Group <- gr[pd$sample, 'Group']
pd$g2 <- paste(pd$drug, pd$Group, sep='_')
pd <- pd[order(nchar(pd$drug)), ]
pd$drug <- factor(pd$drug, levels=unique(pd$drug))
#plot
if(opt$log){
	pd$score=log2(pd$score + 1)
}
p <- ggplot(pd, aes(x=drug, y=score, fill=Group))
if(grepl('violin', opt$type)){
	p <- p + geom_violin()
	if(opt$type=='boxviolin'){
		p <- p + geom_boxplot(aes(group=g2), fill='white', position=position_dodge(0.9), width=0.2, outlier.shape=NA)
	}
	gp.stat <- ggplot_build(p)$data[[1]]
	gp.stat$ymin <- gp.stat$ymax
}else if(opt$type=='boxplot'){
	if(opt$errorbar){
		p <- p + geom_boxplot(linetype="dashed", outlier.shape=NA) + stat_boxplot(aes(ymin=..lower..,ymax=..upper..), size=0.7, outlier.shape=NA)
		p <- p + stat_boxplot(geom = "errorbar",aes(ymin=..ymax..), width=0.4, position=position_dodge(0.75))
		p <- p + stat_boxplot(geom = "errorbar",aes(ymax=..ymin..), width=0.4, position=position_dodge(0.75))
	}else{
		p <- p + geom_boxplot(outlier.size=0.5, size=0.3)
	}
	gp.stat <- ggplot_build(p)$data[[1]]
}
ym <- gp.stat %>% group_by(x=round(x)) %>% summarise(ymax=max(ymax), ymin=min(ymin)) %>% data.frame()
td <- merge(dtr, ym, by='x')
dif <- max(td$ymax) - min(td$ymin)
p <- p + geom_text(data=td, aes(x=drug, y=max(max(pd$score)+dif*0.05), label=lab), inherit.aes=F)
#p <- p + geom_segment(data=td, aes(x=x-0.2, xend=x+0.2, y=ymax+dif*0.01, yend=ymax+dif*0.01), inherit.aes=F)
#p <- p + geom_segment(data=td, aes(x=x-0.2, xend=x-0.2, y=ymax+dif*0.01, yend=ymax), inherit.aes=F)
#p <- p + geom_segment(data=td, aes(x=x+0.2, xend=x+0.2, y=ymax+dif*0.01, yend=ymax), inherit.aes=F)
p <- p + scale_fill_manual(values=colors, guide = guide_legend(title=NULL, nrow=1))
if(opt$switch){
	p <- p + labs(caption=labe.caption)
}
p <- p + coord_cartesian(ylim=c(min(ym$ymin), max(max(pd$score)+dif*0.05))) + labs(x=opt$x, y=opt$y) + theme_classic()
if(opt$angle %in% c(0, 180, 360)){
	xhj <- 0.5
	xvj <- 0.5
}else if(opt$angle == 270){
	xvj <- 0.5
	xhj <- 0
}else if(opt$angle == 90){
	xvj <- 0.5 
	xhj <- 1
}else if(opt$angle < 180){
	xhj <- 1
	xvj <- 1
}else{
	xhj <- 0
	xvj <- 0
}

p <- p + theme(axis.title=element_text(size=20), axis.text=element_text(size=15, color='black'), 
			   axis.text.x=element_text(angle=opt$angle, hjust=xhj, vjust=xvj), legend.position='top', legend.text=element_text(size=12))
ggsave(paste(opt$out, opt$type, '.png', sep=''), plot=p, width=opt$width, height=opt$height, units='in', dpi=600)
ggsave(paste(opt$out, opt$type, '.pdf', sep=''), plot=p, width=opt$width, height=opt$height, units='in')






