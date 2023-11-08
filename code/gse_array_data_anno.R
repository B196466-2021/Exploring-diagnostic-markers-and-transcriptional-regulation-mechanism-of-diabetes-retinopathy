#! /usr/bin/Rscript

##########################
#main programe

file<-strsplit(opt$file, ',')[[1]]
fl<-read.table(file[1], header=T, sep='\t', stringsAsFactors=F, na.strings = c("NA", "null", "", "NULL"), check.names=F)
fl[,1] <- as.character(fl[,1])
if(length(file) > 1){
	for(i in 2:length(file)){
		fcl<-read.table(file[i], header=T, sep='\t', na.strings = c("NA", "null", "", "NULL"))
		fcl[,1] <- as.character(fcl[,1])
		fl<-merge(fl, fcl, by=1)
	}
}
#fma <- as.matrix(fl[,-1])
#fma[is.na(fma)] <- NULL
#fl <- data.frame(ID=fl[,1], fma, stringsAsFactors=F)

anno<-read.table(opt$anno, header=T, sep='\t', comment.char='', quote='', stringsAsFactors=F, check.names=F)
anno[,1] <- as.character(anno[,1])
nsn<-grep('gene symbol|gene_symbol|genesymbol', names(anno), ignore.case=T)
nsgm<-'Gene_Name'

if(length(nsn)==0){
	nsn<-grep('gene name|gene_name|genename', names(anno), ignore.case=T)
	nsgm<-'Gene_Name'
}
if(length(nsn)==0){
	nsn<-grep('miRNA_ID|miRNA ID|miRNA\\t|miRNA$', names(anno), ignore.case=T)
	nsgm<-'miRNA_ID'
}
if(length(nsn)==0){
	nsn<-grep('circRNA_ID|circRNA ID|circRNA', names(anno), ignore.case=T)
	nsgm<-'circRNA_ID'
}
if(length(nsn)==0){
	nsn<-grep('lncRNA', names(anno), ignore.case=T)
	nsgm<-'lncRNA_ID'
}
if(length(nsn)==0){
	nsn<-grep('symbol', names(anno), ignore.case=T)
	nsgm<-'Gene_Name'
}

if(length(nsn)==0){
	nsn<-grep('gene|name', names(anno), ignore.case=T)
	nsgm<-'Gene_Name'
}
rfg<-grep('Regulatory_Feature_Group|Regulatory Feature Group', names(anno), ignore.case=T)
if(length(rfg)!=0){
	nsn<-c(nsn,rfg)
}
ann<-data.frame(ID=anno[,1], anno[,nsn], check.names=F)
if(ncol(ann)==2 && length(rfg)==0){
	names(ann)[2]<-nsgm
}else if(ncol(ann)==3 && length(rfg)!=0){
	names(ann)[c(2,3)]<-c(nsgm, 'Feature_Group')
}
ann<-ann[ann[,2]!='',]
ann<-ann[ann[,2]!='null', ]
ann<-ann[ann[,2]!='NULL', ]
if(opt$ugn){
	gl <- apply(ann[,1:2], 1, function(x){od=data.frame(id=x[1], Gene_Name=strsplit(as.character(x[2]), '\\s*/+\\s*|\\s*;\\s*|\\s*,\\s*|\\|', perl=T)[[1]], stringsAsFactors=F, check.names=F)})
	ann <- do.call(rbind, gl)
}
ann<-ann[ann[,2]!='',]
names(ann)[2]<-nsgm
da<-merge(ann, fl, by=1)
da[is.na(da)]<-0
da <- da[! grepl('^\\d+-|^\\d+$|^\\-+$', da[,2]),]
if(opt$ugn){
	if(any(duplicated(da[, 2]))){
		da <- aggregate(da[, -c(1, 2)], by=list(da[, 2]), mean)
	}else{
		da <- da[, -1]
	}
	names(da)[1] <- nsgm
	if(opt$log){
		da <- data.frame(Gene_Name=da[,1], log2(as.matrix(da[,-1])+1), stringsAsFactors=F, check.names=F)
	}
}
write.table(da, file=opt$out, quote=F, sep='\t', row.names=F)
