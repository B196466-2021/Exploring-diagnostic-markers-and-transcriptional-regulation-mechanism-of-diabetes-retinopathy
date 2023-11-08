#! /usr/bin/Rscript

usage<-function(){
	cat("This script is used for differential expression analysis of array data in GEO\n",
getopt(spec,usage=TRUE),
"Options:
	-d, --data	data matrx (data1,data2,..).
	-p, --phen	phenotype file containing all group informations.
	-c, --comp	compare informations file.
	-o, --out	out directory (default ./).
	-l, --logt	the data have translated by log2 or not (default not translated).
	-t, --type	data type: arry, RNA-seq (default arry).
	-n, --norm	normalize or not (default FALSE).
	-v, --version	display version information.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1)
}

##########################
#get options

library(getopt)
spec = matrix(c(
'data','d',1,'character',
'phen','p',1,'character',
'comp','c',1,'character',
'type','t',1,'character',
'out','o',1,'character',
'logt','l',0,'logical',
'norm','n',2,'logical',
'version','v',0,'logical',
'help','h',0,'logical'
),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$data) || is.null(opt$phen) || is.null(opt$comp) || !is.null(opt$help)) {usage()}

if (is.null(opt$logt)) {opt$logt=F}
if (is.null(opt$out)) {opt$out='./'}
if (is.null(opt$norm)) {opt$norm=F}
if (is.null(opt$type)) {opt$type='arry'}
##########################
#main programe
library(foreach)
library(limma)
options(check.names=F)

pheno<-read.table(opt$phen, header=T, sep='\t', stringsAsFactors=F)
process_pheno <-function(x,pheno){
	        Group<-unlist(strsplit(as.character(pheno[x,2]),","))
        d<-data.frame(Sample=rep(pheno[x,1],length(Group)), Group, stringsAsFactors=F, check.names=F)
		        return(d)
}
ga<-foreach(x=1:nrow(pheno), .combine='rbind') %do% process_pheno(x, pheno)

log2trans=function(ex){
	#qx=as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
	#LogC=(qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	if(opt$logt==F){
		ex[which(ex <= 0)] = NaN
		ex=log2(ex)
		print("data log2 transformed")
	}
	return(ex)
}

files<-strsplit(opt$data, ',')[[1]]
file<-read.table(files[1], header=T, sep='\t', check.names=F, stringsAsFactors=F, comment.char='', quote='', na.strings=c('NA', 'null', 'NULL', 'na'))
file[,1] <- as.character(file[,1])
if(class(file[,2])=="character"){
	gl <- apply(file[,1:2], 1, function(x){od=data.frame(id=x[1], Gene_Name=strsplit(as.character(x[2]), '\\s*/+\\s*|\\s*;\\s*|\\s*,\\s*', perl=T)[[1]], stringsAsFactors=F, check.names=F)})
	gd <- do.call(rbind, gl)
	file2 <- merge(gd, file[, -2], by=1)
	names(file2)[2] <- names(file)[2]
	file <- file2[, -1]
}

if(opt$norm){
	logt <- try(backgroundCorrect(as.matrix(file[, names(file) %in% pheno$Sample]), method="normexp", normexp.method="rma"))
}else{
	logt <- as.matrix(file[, names(file) %in% pheno$Sample])
}
logt2 <- logt[! apply(logt, 1, sum)==0, ]
fl<-data.frame(file[! apply(logt, 1, sum)==0, 1], logt2, check.names=F)
if(any(duplicated(fl[,1]))){
	fl<-aggregate(fl[,-1], by=list(fl[,1]), mean)
}
if(length(files)>1){
	for(i in 2:length(files)){
		file<-read.table(files[i], header=T, sep='\t', row.names=1, comment.char='', check.names=F, stringsAsFactors=F, quote='', na.strings=c('NA', 'null', 'NULL', 'na'))
		if(class(file[,2])=="character"){
			gl <- apply(file[,1:2], 1, function(x){od=data.frame(id=x[1], Gene_Name=strsplit(as.character(x[2]), '\\s*/+\\s*|\\s*;\\s*|\\s*,\\s*', perl=T)[[1]], stringsAsFactors=F, check.names=F)})
			gd <- do.call(rbind, gl)
			file2 <- merge(gd, file[, -2], by=1)
			names(file2)[2] <- names(file)[2]
			file <- file2[, -1]
		}
		if(opt$norm){
			logt <- try(backgroundCorrect(as.matrix(file[, names(file) %in% pheno$Sample]), method="normexp", normexp.method="rma"))
		}else{
			logt <- as.matrix(file[, names(file) %in% pheno$Sample])
		}
		logt2 <- logt[! apply(logt, 1, sum)==0, ]
		fls <- data.frame(file[! apply(logt, 1, sum)==0, 1], logt2, check.names=F)
		if(any(duplicated(fls[,1]))){
			fls<-aggregate(fls[,-1], by=list(fls[,1]), mean)
		}
		fl<-merge(fl, fls, by=1)

	}
}
rownames(fl)<-fl[,1]
fl<-fl[,-1]
if(opt$type=='arry'){
	boxda<-data.frame(expression=as.numeric(as.matrix(fl)), sample=rep(colnames(fl), each=nrow(fl)), check.names=F)
	pdf(paste(opt$out, 'pre_normalize_bosplot.pdf', sep=''))
	par(las=2)
	boxplot(expression ~ sample, data=boxda,  col=rainbow(length(colnames(fl))), outline=F, xlab=NULL)
	dev.off()
	fl <- normalizeBetweenArrays(as.matrix(fl), method='quantile')
	boxda<-data.frame(expression=as.numeric(as.matrix(fl)), sample=rep(colnames(fl), each=nrow(fl)), check.names=F)
	pdf(paste(opt$out, 'pos_normalize_bosplot.pdf', sep=''))
	par(las=2)
	boxplot(expression ~ sample, data=boxda,  col=rainbow(length(colnames(fl))), outline=F, xlab=NULL)
	dev.off()
}
allexp<-data.frame(rownames(fl), fl, check.names=F)
names(allexp)[1]<-names(file)[1]
write.table(allexp, file=paste(opt$out, '/gene_expression.txt', sep=''), quote=F, sep='\t', row.names=F)
comp<-read.table(opt$comp,header=T,sep="\t",colClasses = c(rep("character", 2), rep("numeric", 3)))

if(! dir.exists(paste(opt$out,'/sep',sep=''))){
	dir.create(paste(opt$out,'/sep',sep=''), recursive=T)
}
for(i in 1:nrow(comp)){
	if(! (comp$Test[i] %in% ga$Group && comp$Control[i] %in% ga$Group)){
		next
	}
	sug<-ga[ga$Group %in% c(comp$Test[i], comp$Control[i]),]
	rownames(sug)<-sug$Sample
	expt<-data.frame(t(fl[, colnames(fl) %in% sug$Sample]), check.names=F)
	gexp<-merge(sug, expt, by=0)
	rownames(gexp)<-gexp[,1]
	expm<-data.frame(t(gexp[,-c(1,2,3)]), check.names=F)
	expm <- expm[! apply(expm, 1, sum)==0, ]
	expm <- data.frame(apply(expm, 2, function(x){x[x==0]= min(x[!x==0]); return(x)}), check.names=F)
	lot<-log2trans(expm)
	expm<-lot
	sug<-gexp[,c(2,3)]
	design<-model.matrix(~ Group + 0, sug)
	colnames(design)<-levels(factor(sug$Group))
	rownames(design)<-sug$Sample
	fit<-lmFit(expm, design)
	cont.matrix <- makeContrasts(contrasts=paste(comp$Test[i], '-', comp$Control[i], sep=''), levels=design)
	fit2<-contrasts.fit(fit, cont.matrix)
	if(all(table(sug$Group)<2)){
		fiteBayes<-data.frame(fit2$coefficients, P.Value=1, adj.P.Val=1, check.names=F)
		names(fiteBayes)[1] <- 'logFC'
	}else{
		fit3<-eBayes(fit2)
		fiteBayes<-topTable(fit3, coef=1, number=100000000, sort.by="logFC")
	}
	out_names<-c('logFC', 'P.Value', 'adj.P.Val')
	oall<-merge(fiteBayes[,colnames(fiteBayes)%in% out_names], fit$coefficients, by=0)
	rownames(oall)<-oall[,1]
	if(!all(colnames(expm) %in% colnames(oall))){
		oall<-merge(oall[,-1], expm, by=0)
	}
	oall<-data.frame(oall[,1:2], FC=2**oall$logFC, oall[,-(1:2)], check.names=F)
	colnames(oall)[1:5]<-c('Gene_Name', 'log2FC', 'FC',  'p_value', 'q_value')
	#if(logtr==F && opt$logt==F){
	#	oall$log2FC<-log2(oall[, comp$Test[i]]/oall[, comp$Control[i]])
	#}
	upBayes<-oall[oall$log2FC>=log2(comp$fc[i]) & oall$p_value <=comp$p_value[i] & oall$q_value<=comp$q_value[i],]
	upBayes <- upBayes[order(upBayes$log2FC, decreasing=T),]
	if(nrow(upBayes)>0){
		upg<-data.frame(unique(upBayes[,1]), check.names=F)
		write.table(upg, file=paste(opt$out, '/sep/up_', comp$Test[i], '_vs_', comp$Control[i], '_gene.txt', sep=''), quote=F, sep='\t', row.names=F, col.names=F)
		write.table(upBayes, file=paste(opt$out, '/sep/up_', comp$Test[i], '_vs_', comp$Control[i], '.txt', sep=''), quote=F, sep='\t', row.names=F)
	}else{
		print('Warring: no up regulated data.')
	}
	downBayes<-oall[oall$log2FC<=log2(1/comp$fc[i]) & oall$p_value <=comp$p_value[i] & oall$q_value<=comp$q_value[i],]
	downBayes <- downBayes[order(downBayes$log2FC),]
	if(nrow(downBayes)>0){
		downg<-data.frame(unique(downBayes[,1]), check.names=F)
		write.table(downg, file=paste(opt$out, '/sep/down_', comp$Test[i], '_vs_', comp$Control[i], '_gene.txt', sep=''), quote=F, sep='\t', row.names=F, col.names=F)
		write.table(downBayes, file=paste(opt$out, '/sep/down_', comp$Test[i], '_vs_', comp$Control[i], '.txt', sep=''), quote=F, sep='\t', row.names=F)
	}else{
		print('Warring: no down regulated data.')
	}
	write.table(oall, file=paste(opt$out, '/', comp$Test[i], '_vs_', comp$Control[i], '.txt', sep=''), quote=F, sep='\t', row.names=F)
	if(nrow(upBayes)>0 && nrow(downBayes)>0){
		updown <- rbind(upBayes, downBayes[order(downBayes$log2FC, decreasing=T),])
		write.table(updown, file=paste(opt$out, '/sep/up+down_', comp$Test[i], '_vs_', comp$Control[i], '.txt', sep=''), quote=F, sep='\t', row.names=F)
	}
}
	
