#! /usr/bin/Rscript

##########################
#main programe
library(VennDiagram)
library(gplots)
library(venn)
library(ggplot2)

files<-strsplit(opt$file, ',')[[1]]

if(is.null(opt$pos) && length(files)==3){
	pos<-c(0, 0, 180)
}else if(is.null(opt$pos) && length(files)==2){
	pos<-c(0, 0)
}else if(!is.null(opt$pos)){
	pos<-as.numeric(strsplit(opt$pos, ',')[[1]])
}

fl<-list()
if(length(files)>1){
	if(is.null(opt$name)){
		stop('\033[31mYou must provide the name for each list file.\033[0m')
	}
	ns<-strsplit(opt$name, ',')[[1]]
	if(length(files) != length(ns)){
		stop('\033[31mThe number of names not equal to list files.\033[0m')
	}
	for(i in 1:length(files)){
		fl[[i]]<-unique(readLines(files[i]))
	}
	names(fl)<-ns
}else{
	fd<-read.table(files[1], header=T, row.names=1, sep='\t', stringsAsFactors=F)
	for (j in 1:ncol(fd)){
		ronm<-rownames(fd)
		fl[[j]]<-unique(ronm[fd[, j]>0])
	}
	names(fl)<-names(fd)
}
insr<-venn(fl,show.plot=F,intersections=T)
iar<-attr(insr,"intersections")
nsp<-strsplit(names(iar),':')
lnm<-0
for(i in 1:length(nsp)){
	if(length(nsp[[i]])==length(fl)){
		lnm<-i
	}
}
if(is.null(opt$color)){
	colors <- rainbow(length(fl))
}else{
	colors <- strsplit(opt$color, ',')[[1]]
}
if(lnm>0){
	ino<-data.frame(igen=iar[[lnm]])
	write.table(ino, file=paste(opt$out, 'gene.txt', sep=''), quote=F, sep='\t', row.names=F, col.names=F)
}
for(i in 1:length(fl)){
	ui <- fl[[i]][!fl[[i]] %in% unlist(fl[-i])]
	writeLines(ui, paste(names(fl[i]), '_only.txt', sep=''))
}
if(length(fl)==2){
	p <- venn.diagram(fl, filename=NULL, scaled=F, cat.pos=pos, cat.dist=0.03, fill=colors, col=NA, cat.cex=1.5, cex=2, alpha=0.3)
	#uf1 <- fl[[1]][!fl[[1]] %in% fl[[2]]]
	#uf2 <- fl[[2]][!fl[[2]] %in% fl[[1]]]
	#writeLines(uf1, paste(names(fl[1]), '_only.txt', sep=''))
	#writeLines(uf2, paste(names(fl[2]), '_only.txt', sep=''))
}else if(length(fl) > 5){
	p <- venn::venn(fl, zcolor ='style', ggplot=T, col='gray60', zcolor=paste(colors, collapse=', '))
}else if(! is.null(opt$pos)){
	p <- venn.diagram(fl, filename=NULL, scaled=F, fill=colors, alpha=0.3, cat.pos=pos, col=NA, cat.cex=1.5, cex=1.5, imagetype='png',margin=c(0.1,0.1,0.1,0.1))
}else{
	p <- venn.diagram(fl, filename=NULL, scaled=F, fill=colors, alpha=0.3, col=NA, cat.cex=1.5, cex=1.5, imagetype='png',margin=c(0.1,0.1,0.1,0.1))
}
if(class(p)=='gList'){
	png(paste(opt$out, 'venn_plot.png', sep=''), width=5, height=5, units='in', res=600)
	grid.draw(p)
	dev.off()
	pdf(paste(opt$out, 'venn_plot.pdf', sep=''), width=5, height=5)
	grid.draw(p)
	dev.off()
}else{
	ggsave(paste(opt$out, 'venn_plot.png', sep=''), plot=p, width=5, height=5, units='in', dpi=600)
	ggsave(paste(opt$out, 'venn_plot.pdf', sep=''), plot=p, width=5, height=5, units='in')
}
file.remove(list.files(path='./', pattern='VennDiagram.*.log|Rplots.*.pdf'))

