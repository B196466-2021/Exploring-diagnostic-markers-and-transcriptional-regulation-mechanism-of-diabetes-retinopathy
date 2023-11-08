#! /usr/bin/Rscript

##########################
#main programe
library(psych)
library(ggplot2)
library(parallel)
f1 <- read.table(opt$f1, header=T, sep='\t', stringsAsFactors=F, check.names=F)
if(!class(f1[, 1]) %in% c('numeric', 'integer')){
	rownm <- f1[, 1]
	conm1 <- names(f1)[-1]
	f1 <- data.frame(f1[, -1])
	rownames(f1) <- rownm
	names(f1) <- conm1
}else{
	 conm1 <- names(f1)
}
if(is.null(opt$f2)){
	cor.re <- corr.test(f1, adjust='none')
}else{
	f2 <- read.table(opt$f2, header=T, sep='\t', stringsAsFactors=F, check.names=F)
	if(!class(f2[, 1]) %in% c('numeric', 'integer')){
		rownm <- f2[, 1]
		conm2 <- names(f2)[-1]
		f2 <- data.frame(f2[, -1])
		rownames(f2) <- rownm
		names(f2) <- conm2 
		if(length(which(rownames(f2) %in% rownames(f1))) >2){
			rnams <- rownames(f2)[rownames(f2) %in% rownames(f1)]
			f2 <- data.frame(f2[rnams, ])
			names(f2) <- conm2
			f1 <- data.frame(f1[rnams, ])
			rownames(f2) <- rownames(f1) <- rnams
			names(f1) <- conm1
		}else if(length(which(colnames(f1)  %in%  colnames(f2))) >2){
			f1 <- t(f1)
			f2 <- t(f2)
			rownm <- rownames(f2)[rownames(f2) %in% rownames(f1)]
			conm1 <- colnames(f1)
			conm2 <- colnames(f2)
			f1 <- data.frame(f1[rownm, ])
			colnames(f1) <- conm1
			f2 <- data.frame(f2[rownm, ])
			colnames(f2) <- conm2
			rownames(f1) <- rownames(f2) <- rownm
			warnings('translate row to column of f1 and f2')
		}else if(length(which(colnames(f1)  %in%  rownames(f2))) >2){
			f1 <- t(f1)
			rownm <- rownames(f2)[rownames(f2) %in% rownames(f1)]
			conm1 <- colnames(f1)
			conm2 <- colnames(f2)
			f1 <- data.frame(f1[rownm, ])
			colnames(f1) <- conm1
			f2 <- data.frame(f2[rownm, ])
			colnames(f2) <- conm2
			rownames(f1) <- rownames(f2) <- rownm
			warnings('translate row to column of f1')
		}else if(length(which(colnames(f2)  %in%  rownames(f1))) >2){
			f2 <- t(f2)
			rownm <- rownames(f2)[rownames(f2) %in% rownames(f1)]
			conm1 <- colnames(f1)
			conm2 <- colnames(f2)
			f1 <- data.frame(f1[rownm, ])
			colnames(f1) <- conm1
			f2 <- data.frame(f2[rownm, ])
			colnames(f2) <- conm2
			rownames(f1) <- rownames(f2) <- rownm
			warnings('translate row to column of f2')
		}else{
			stop('row names of tatle 1 are not same as tatle 2.')
		}
	}else{
		conm2 <- names(f2)
	}
	if(all(c(ncol(f1), ncol(f2))==1)){
		data <- cbind(f1, f2)
		p <- corr_scatter(x=conm1, y=conm2, data=data)
		ggsave(paste(opt$out, 'scatterhist_plot.pdf', sep=''), plot=p, width=opt$width, height=opt$height, units='in')
		ggsave(paste(opt$out, 'scatterhist_plot.png', sep=''), plot=p, width=opt$width, height=opt$height, units='in', dpi=600)
		file.remove('Rplots.pdf')
		quit(save='no')
	}else if(any(c(ncol(f1), ncol(f2))==1)){
		if(ncol(f1)==1){
			x=f2
			y=f1
		}else{
			x=f1
			y=f2
		}
		coret <- corr_poiont_plot(x, y)
		ggsave(paste(opt$out, 'correlation_plot.pdf', sep=''), plot=coret$plot, width=opt$width, height=opt$height, units='in')
		ggsave(paste(opt$out, 'correlation_plot.png', sep=''), plot=coret$plot, width=opt$width, height=opt$height, units='in', dpi=600)
		write.table(coret$cor, paste(opt$out, 'correlation_result.txt', sep=''), sep='\t', quote=F, row.names=F)
		q(save='no')
	}else{
		mc <-getOption("mc.cores", opt$thread)
		corl <- mclapply(colnames(f2), mcot, da1=f1, da2=f2, mc.cores = mc)
		pl <- rl <- list()
		for(i in 1:length(corl)){
			pl[[i]] <- corl[[i]]$p
			rl[[i]] <- corl[[i]]$r
		}
		cor.re <- list()
		cor.re[['r']] <- do.call(cbind, rl)
		cor.re[['p']] <- do.call(cbind, pl)
	}
}
rd <- data.frame(ID=rownames(cor.re$r), cor.re$r, check.names=F)
pd <- data.frame(ID=rownames(cor.re$p), cor.re$p, check.names=F)
write.table(rd, file=paste(opt$out, 'correlation_coefficient.txt', sep=''), sep='\t', row.names=F, quote=F)
write.table(pd, file=paste(opt$out, 'correlation_p_value.txt', sep=''), sep='\t', row.names=F, quote=F)

if(is.null(opt$f2)){
	cor.pd <- data.frame()
	for(i in 1:(ncol(cor.re$r)-1)){
		sub.pd <- data.frame(x=colnames(cor.re$r)[i], y=rownames(cor.re$r)[-(1:i)], r=cor.re$r[-(1:i), i], p=round(cor.re$p[-(1:i), i], 2), stringsAsFactors=F, check.names=F)
	    cor.pd <- rbind(cor.pd, sub.pd)
	}
}else{
	cor.pd <- data.frame(x=rep(rownames(cor.re$r), ncol(cor.re$r)), y=rep(colnames(cor.re$r), each=nrow(cor.re$r)), r=as.numeric(cor.re$r), p=round(as.numeric(cor.re$p), 2), stringsAsFactors=F)
}
#if(max(nchar(rownames(cor.re$r))) > max(nchar(colnames(cor.re$r)))){
#	names(cor.pd)[1:2] <- c('y', 'x')
#	cor.pd$x <- factor(cor.pd$x, levels=colnames(cor.re$r))
#	cor.pd$y <- factor(cor.pd$y, levels=rev(rownames(cor.re$r)))
#}else{
	cor.pd$x <- factor(cor.pd$x, levels=rownames(cor.re$r))
	cor.pd$y <- factor(cor.pd$y, levels=colnames(cor.re$r))
#}
cor.pd$lb <- as.character(cor.pd$p)
cor.pd$lb[is.na(cor.pd$p)|cor.pd$p >=0.05] <- ''
cor.pd$lb[cor.pd$p >=0.01 &cor.pd$p <0.05] <- '*'
cor.pd$lb[cor.pd$p >=0.001 &cor.pd$p <0.01] <- '**'
cor.pd$lb[cor.pd$p <0.001] <- '***'
labe.caption <- expression(paste('*: ', 0.01 <= p, '< 0.05;  **: ', 0.001 <= p, '< 0.01;  ***: p < 0.001'))
na.r <- cor.pd$r[!is.na(cor.pd$r)]
colors <- strsplit(opt$color, ',')[[1]]
if(min(na.r) >= 0){
	colors <- colors[2:3]
	values <- scales::rescale(c(min(na.r), max(na.r)))
}else if(max(na.r) <=0){
	colors <- colors[1:2]
	values <- scales::rescale(c(min(na.r), max(na.r)))
}else{
	colors <- colors
	values <- scales::rescale(c(min(na.r), 0,  max(na.r)))
}

p <- ggplot(cor.pd, aes(x, y, fill=r)) + geom_tile(colour='grey50')
p <- p + scale_fill_gradientn(colors=colors, values=values, na.value='gray70', guide=guide_colourbar(title='R', title.position='top', title.hjust=0.2, ticks=F))
p <- p + scale_y_discrete(expand=c(0, 0)) + scale_x_discrete(, expand=c(0, 0))
if(opt$pval){
	p <- p + geom_text(aes(label=lb), size=5, vjust=0.8)
	p <- p + labs(caption=labe.caption)
}
p <- p + theme(axis.title=element_blank(), axis.ticks.length=unit(0, 'mm'), legend.background=element_blank(), panel.grid=element_blank(), panel.background=element_blank(), panel.border = element_rect(fill = NA,colour = "grey20"))
if(opt$name){
	p <- p + theme(axis.text=element_text(size=12, colour='black'), axis.text.x=element_text(angle=45, hjust=1))
}else{
	p <- p + theme(axis.text=element_blank())
}
if(is.null(opt$f2)){
	p <- p + theme(legend.position=c(0.8, 0.2))
}else{
p <- p + theme(legend.position='right')
}   

ggsave(paste(opt$out, 'correlation_heatmap.pdf', sep=''), plot=p,  width=opt$width, height=opt$height, units='in')
ggsave(paste(opt$out, 'correlation_heatmap.png', sep=''), plot=p,  width=opt$width, height=opt$height, units='in', dpi=600)

names(cor.pd) <- c('y', 'x', 'R', 'p_value')
write.table(cor.pd[,1:4], file=paste(opt$out, 'correlation_result.txt', sep=''), sep='\t', quote=F, row.names=F)
