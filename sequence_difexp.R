#! /usr/bin/Rscript
usage<-function(){
	cat("This script is used for differential expression analysis of circleRNA.\n",
getopt(spec,usage=TRUE),
"Options:
	-f, --file	file containing count of all samples.
	-p, --pheno	group information file.
	-c, --comp	compare information file.
	-t, --type	the expression data type (counts, FPKM, RPKM)
	-m, --cpm	the minimum value of cpm for filter(default 10).
	-o, --out	output directory(default ./).
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
'pheno','p',1,'character',
'comp','c',1,'character',
'type','t',1,'character',
'circRNA','r',0,'logical',
'cpm','m',1,'double',
'out','o',1,'character',
'version','v',0,'logical',
'help','h',0,'logical'
),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$file) || is.null(opt$pheno) || is.null(opt$comp) || !is.null(opt$help)) {usage()}

if (is.null(opt$type)) {stop('You must provide the data type by -t or --type.')}
if (is.null(opt$cpm)) {opt$cpm=10}
if (is.null(opt$out)) {opt$out="./"}

if(! dir.exists(opt$out)){
	dir.create(opt$out, recursive=T)
}
##########################
#main programe

library(statmod)
library(methods)
library(locfit)
library(foreach)

group<-read.table(opt$pheno,header=T,sep="\t",colClasses="character")
comp<-read.table(opt$comp,header=T,sep="\t",colClasses = c(rep("character", 2), rep("numeric", 3)))
if(! is.null(opt$anno)){
	anno<-read.table(opt$anno,header=T,sep="\t",check.names =F)
}

process_pheno <-function(x,pheno){
	Group<-unlist(strsplit(as.character(pheno[x,2]),","))
	d<-data.frame(Sample=rep(pheno[x,1],length(Group)),Group, stringsAsFactors=F)
	return(d)
}
group<-foreach(x=1:nrow(group),.combine='rbind') %do% process_pheno(x,group)

rawdata <- read.delim(opt$file, check.names=FALSE, stringsAsFactors=FALSE)
group <- group[group$Sample %in% names(rawdata),]
rawdata <- rawdata[, c(names(rawdata)[1], group$Sample)]
if(any(duplicated(rawdata[, 1]))){
	mdata<-aggregate(rawdata[, -1], list(rawdata[,1]), mean)
}else{
	mdata <- rawdata
}
rawdata <- mdata[, -1]
rownames(rawdata) <- mdata[, 1]
sr<-apply(rawdata,1,sum)

rawdata<-rawdata[sr>0, ]
ga<-group[group$Sample %in% names(rawdata),]
#write.table(group,file=paste(opt$out,"/pheno2.txt",sep=""),quote=F,sep="\t",row.names=F,na ="NA")
#ga<-foreach(x=1:nrow(group),.combine='rbind') %do% process_pheno(x,group)
if(opt$type=='counts'){
	library(edgeR)
	source("/home/dev/DEV-wangdh/exactTest.R")
	comp2<-data.frame()
	for(i in 1:nrow(comp)){
		if(comp$Test[i] %in% ga$Group && comp$Control[i] %in% ga$Group){
			comp2<-rbind(comp2,comp[i,])
		}
	}
	fc<-c()
	for(i in unique(ga$Group)){
		subc<-data.frame(rawdata[,colnames(rawdata) %in% ga[ga$Group==i,]$Sample],check.names=F)
		rownames(subc)<-rownames(rawdata)
		subfc<-data.frame(subc[apply(cpm(subc),1,mean)>opt$cpm,],check.names=F)
		if(nrow(subfc)>0){
			rownames(subfc)<-rownames(subc)[apply(cpm(subc),1,mean)>opt$cpm]
			fc<-c(fc,rownames(subfc))
		}
	}
	
	filterc<-data.frame(rawdata[rownames(rawdata) %in% unique(fc),],check.names=F)
	filterc2<-data.frame(rownames(filterc),filterc,check.names=F)
	names(filterc2)[1]<-"ID"
	write.table(filterc2,file=paste(opt$out,"/gene_expression.txt",sep=""),quote=F,sep="\t",row.names=F,na ="NA")
	deal<-data.frame()
	if(!dir.exists(paste(opt$out,"/sep",sep=""))){
		dir.create(paste(opt$out,"/sep",sep=""))
	}
	rcom<-c()
	for(i in 1:nrow(comp2)){
		if(! all(c(comp$Test[i], comp$Control[i]) %in% ga$Group)){
			warnings(paste('groups ', comp$Test[i], 'and', comp$Control[i], ' are not all in ', opt$pheno, sep=''))
			next
		}
		pair<-c(as.character(comp2$Control[i]),as.character(comp2$Test[i]))
		gs<-ga[which(ga$Group %in% pair),]
		ds<-filterc[,which(colnames(filterc) %in% gs$Sample)]
		sds<-apply(ds,2,sum)
		if(0 %in% sds){
			sds[sds>0]<-5
			sds[sds==0]<-1
			sds[sds==5]<-1
			ds<-rbind(ds,sds)
			rownames(ds)[nrow(ds)]<-'dummy'
			ds<-rbind(ds,sds)
			rownames(ds)[nrow(ds)]<-'dummy2'
		}
		gso<-merge(data.frame(Sample=colnames(filterc)),gs,by=1,sort=F)
		y<-DGEList(counts=ds,group=gso$Group)
		gcpm<-merge(gso,t(cpm(y)),by.x=1,by.y="row.names")
		if(any(duplicated(gcpm$Group))){
			mcpm<-aggregate(gcpm[,3:ncol(gcpm)],list(gcpm$Group),mean)
		}else{
			mcpm<-data.frame(Group=gcpm$Group, gcpm[,3:ncol(gcpm)], stringsAsFactors=F)
		}
		mcpm2<-merge(data.frame(c(pair[2],pair[1]),c("Test_value","Control_value")),mcpm,by=1)
		rownames(mcpm2)<-mcpm2[,2]
		mcpm3<-t(mcpm2[,-c(1,2)])
		
		yf <- y[apply(mcpm3>opt$cpm,1,any),]
		yf$samples$lib.size <- colSums(yf$counts)
	
		#sampleNum <- length(y$samples$group) / length(unique(y$samples$group))
		#keep <- rowSums(cpm(y) > opt$cpm) >= sampleNum
		#yf <- y[keep,]
		#yf$samples$lib.size <- colSums(yf$counts)
		yn <- calcNormFactors(yf)
		gs<-table(yn$samples$group)
		if(length(gs)==1){
			next
		}
		rcom<-c(rcom,i)
		if(1 %in% gs){
			comp2$p_value[i]<-1
			comp2$q_value[i]<-1
			#yd<-estimateGLMCommonDisp(yn,method="deviance", robust=TRUE , subset=NULL)
			dexp<-exactTest(yn,pair=pair,dispersion=0.2^2)
			dex<-topTags(dexp,n=nrow(dexp$table))$table
			dex[,c(3,6)]<-1
			yc=yn
	
		}else{
			#pdf(file = paste(opt$out,"/",pair[2],"_vs_",pair[1],"_MDS.pdf",sep=""))
			#plotMDS(yn)
			#dev.off()
			#td<-yn$counts[c(1:5,nrow(yn$counts)),]
			#print(td)
			#estimateDisp(td, robust = TRUE)
			#print('yes')
			yd <- estimateDisp(yn, robust = TRUE,min.row.sum=0)
			#pdf(file = paste(opt$out,"/",pair[2],"_vs_",pair[1],"_BCV.pdf",sep=""))
			#plotBCV(yd)
			#dev.off()
			dexp<-exactTest(yd,pair=pair)
			dex<-topTags(dexp,n=nrow(dexp$table))$table
			yc=yd
		}
		dt<-data.frame(dex,Status=rep("OK",nrow(dex)),check.names=FALSE)
		m<-merge(yc$counts,dt,by="row.names",all=T)
		s<-as.character(m$Status)
		s[is.na(s)]<-"NOTEST"
		m$Status<-s
		#m2<-data.frame(m,Test=rep(pair[2],nrow(m)),Control=rep(pair[1],nrow(m)))
		m2<-data.frame(m)
		m3<-data.frame(ID=m2[,1],Test=pair[2],Control=pair[1],Status=m2$Status,Test_value=m2$Test_value,Control_value=m2$Control_value,log2FC=m2$logFC,p_value=m2$PValue,q_value=m2$FDR,check.names=FALSE)
		if(is.null(opt$anno)){
			m4<-m3
		}else{
			m4<-merge(anno,m3,by=1)
		}
	
		mde<-merge(dex,cpm(yc),by="row.names")
		mde2<-data.frame(Gene_Name=mde[,1],log2FC=mde$logFC,FC=2**mde$logFC,p_value=mde$PValue,q_value=mde$FDR,mde$Test_value,mde$Control_value,mde[,8:ncol(mde)],check.names=FALSE)
		names(mde2)[c(6,7)]<-c(pair[2],pair[1])
		if(is.null(opt$anno)){
			mde3<-mde2
		}else{
			mde3<-merge(anno,mde2,by=1)
		}
		if (!is.null(mde3$Locus)){
			mde3<-mde3[order(mde3$Locus),]
		}
		write.table(mde3,file=paste(opt$out,"/",pair[2],"_vs_",pair[1],".txt",sep=""),quote=F,sep="\t",row.names=F)
	
		oupde <- subset(mde3, log2FC >= log2(comp2[i, 'fc']) & p_value <= comp2[i, 'p_value'] & q_value <= comp2[i, 'q_value'])
		oupde <- oupde[order(oupde$log2FC, decreasing=T), ]
		odownde <- subset(mde3, log2FC <= -log2(comp2[i, 'fc']) & p_value <= comp2[i, 'p_value'] & q_value <= comp2[i, 'q_value'])
		odownde <- odownde[order(odownde$log2FC), ]
		updown <- rbind(oupde, odownde[order(odownde$log2FC, decreasing=T), ])
		write.table(oupde,file=paste(opt$out,"/sep/up_",pair[2],"_vs_",pair[1],".txt",sep=""),quote=F,sep="\t",row.names=F)
		write.table(odownde,file=paste(opt$out,"/sep/down_",pair[2],"_vs_",pair[1],".txt",sep=""),quote=F,sep="\t",row.names=F)
		write.table(updown,file=paste(opt$out,"/sep/up+down_",pair[2],"_vs_",pair[1],".txt",sep=""),quote=F,sep="\t",row.names=F)
		writeLines(oupde$Gene_Name, con=paste(opt$out,"/sep/up_",pair[2],"_vs_",pair[1],"_gene.txt",sep=""))
		writeLines(odownde$Gene_Name, con=paste(opt$out,"/sep/down_",pair[2],"_vs_",pair[1],"_gene.txt",sep=""))
		
		if(nrow(deal)==0){
			deal<-m4
		}else{
			deal<-rbind(deal,m4)
		}
	}
	comp2<-comp2[rcom,]
	write.table(comp2,file=paste(opt$out,"/compare2.txt",sep=""),quote=F,sep="\t",row.names=F)
	#write.table(deal,file=paste(opt$out,"/all.diff",sep=""),quote=F,sep="\t",row.names=F)
}else if(opt$type %in% c('FPKM', 'RPKM')){
	library(ballgown)
	names(mdata)[1] <- 'Gene_Name'
	write.table(mdata, file=paste(opt$out, "/gene_expression.txt", sep=""),quote=F,sep="\t",row.names=F)
	for(m in 1:nrow(comp)){
		if(! all(c(comp$Test[m], comp$Control[m]) %in% ga$Group)){
			warning(paste('groups ', comp$Test[m], 'and', comp$Control[m], ' are not all in ', opt$pheno, sep=''))
			next
		}
		subg <- ga[ga$Group %in% c(comp[m, 'Test'], comp[m, 'Control']),]
		subg$Group[subg$Group==comp[m, 'Test']] <- 1
		subg$Group[subg$Group==comp[m, 'Control']] <- 0
		subexp <- rawdata[, subg$Sample]
		deg <- stattest(gowntable=subexp, pData=subg, meas='FPKM', covariate = 'Group', getFC = T, log = T, libadjust = FALSE, feature='gene')
		mexp <- as.matrix(subexp)
		log2exp <- log2(mexp+1)
		mean.test <- apply(log2exp[, subg$Sample[subg$Group==1]], 1, mean)
		mean.control <- apply(log2exp[, subg$Sample[subg$Group==0]], 1, mean)
		mean.tc <- data.frame(mean.test, mean.control)
		names(mean.tc) <- c(comp[m, 'Test'], comp[m, 'Control'])
		if(!dir.exists(paste(opt$out,"/sep",sep=""))){
			dir.create(paste(opt$out,"/sep",sep=""))
		}
		deg.all <- data.frame(Gene_Name=deg$id, log2FC=log2(deg$fc), FC=deg$fc, p_value=deg$pval, q_value=deg$qval, mean.tc, log2exp, check.names=F)
		write.table(deg.all, file=paste(opt$out, "/", comp[m, 'Test'], "_vs_", comp[m, 'Control'], ".txt", sep=""), quote=F, sep="\t", row.names=F)
		oupde <- subset(deg.all, log2FC >= log2(comp[m, 'fc']) & p_value <= comp[m, 'p_value'] & q_value <= comp[m, 'q_value'])
		oupde <- oupde[order(oupde$log2FC, decreasing=T), ]
		write.table(oupde, file=paste(opt$out, "/sep/up_", comp[m, 'Test'] ,"_vs_", comp[m, 'Control'], ".txt", sep=""),quote=F,sep="\t",row.names=F)
		writeLines(as.character(oupde$Gene_Name), con=paste(opt$out, "/sep/up_", comp[m, 'Test'], "_vs_", comp[m, 'Control'], "_gene.txt", sep=""))
		odownde <- subset(deg.all, log2FC <= -log2(comp[m, 'fc']) & p_value <= comp[m, 'p_value'] & q_value <= comp[m, 'q_value'])
		odownde <- odownde[order(odownde$log2FC), ]
		write.table(odownde,file=paste(opt$out, "/sep/down_", comp[m, 'Test'], "_vs_", comp[m, 'Control'], ".txt", sep=""), quote=F, sep="\t", row.names=F)
		writeLines(as.character(odownde$Gene_Name), con=paste(opt$out, "/sep/down_", comp[m, 'Test'], "_vs_", comp[m, 'Control'], "_gene.txt", sep=""))
		updown <- rbind(oupde, odownde[order(odownde$log2FC, decreasing=T), ])
		write.table(updown,file=paste(opt$out, "/sep/up+down_", comp[m, 'Test'], "_vs_", comp[m, 'Control'], ".txt", sep=""), quote=F, sep="\t", row.names=F)
	}
}
