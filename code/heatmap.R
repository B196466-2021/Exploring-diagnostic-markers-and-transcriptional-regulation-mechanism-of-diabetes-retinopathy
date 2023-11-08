#! /usr/bin/Rscript

##########################
#main programe

if (dir.exists(opt$dir) || (! is.null(opt$comp) && file.exists(opt$dir))){
	if(is.null(opt$comp)) {
		stop("Can not find compare information file.")
	}
	comp<-read.table(opt$comp,header=T,sep="\t",check.names=F,colClasses = c(rep("character", 2), rep("numeric", 3)))

	for(i in 1:nrow(comp)){
		if(! all(c(comp$Test[i], comp$Control[i]) %in% ga$Group)){
			warnings(paste('groups ', comp$Test[i], 'and', comp$Control[i], ' are not all in ', opt$pheno, sep=''))
			next
		}
		values<-c("p_value","q_value")
		value<-min(c(comp$p_value[i],comp$q_value[i]))
		if(comp$p_value[i]<comp$q_value[i]){
			value<-comp$p_value[i]
			pq<-"p_value"
		}else{
			pq<-"q_value"
			value<-comp$q_value[i]
		}
		if(dir.exists(opt$dir) && file.exists(opt$comp)){
			file<-paste(opt$dir,"/",comp$Test[i],"_vs_",comp$Control[i],".txt",sep="")
		}else if(file.exists(opt$dir)){
			file<-opt$dir
		}
		if(!file.exists(file)){
			next
		}
		dife<-read.table(file,header=T,sep="\t",row.names=opt$name, check.names=F, stringsAsFactors=F, quote="")
		if(!is.null(opt$type)){
			lncRNA<-c("processed_transcript","lincRNA","3prime_overlapping_ncrna","antisense","non_coding","sense_intronic","sense_overlapping","TEC","known_ncrna","macro_lncRNA","bidirectional_promoter_lncrna")
			sncRNA<-c("snRNA","snoRNA","rRNA","Mt_tRNA","Mt_rRNA","misc_RNA","miRNA","ribozyme","sRNA","scaRNA","vaultRNA")
		#	dife<-read.table(file,header=T,sep="\t",row.names=1,check.names=F)
			otp<-c("Gene_Type","Trans_Type")
			tm<-otp[which(otp %in% names(dife))]
			if(length(tm)==0){
				tm<-"Gene_Type"
			}
			if(opt$type == "all"){
				type<-type<-unique(dife[,eval(tm)])
			}else{
				type<-unlist(strsplit(opt$type,","))
				if("lncRNA" %in% type){
					type<-c(type,lncRNA)
				}
				if("sncRNA" %in% type){
					type<-c(type,sncRNA)
				}
			}
	
			if(any(type %in% dife[,eval(tm)])){
				dif<-dife[which(dife[,eval(tm)] %in% type),]
			}else{
				stop("The type provided doesn't exist.\n")
			}
		}else{
			dif<-dife
		}
		vs<-as.character(t(comp[i,1:2]))
		g<-ga[which(ga$Group %in% vs),]
		g$Group <- factor(as.character(g$Group), levels=as.character(c(comp$Test[i], comp$Control[i])))
		g <- g[order(as.numeric(g$Group)),]
		anno_col<-data.frame(Group=g[,2])
		anno_col$Group <- factor(as.character(anno_col$Group), levels=as.character(c(comp$Test[i], comp$Control[i])))
		anno_colors <- list(Group=get(paste('pal_', opt$palette, sep=''))()(length(unique(anno_col$Group)))[as.numeric(unique(anno_col$Group))])
		names(anno_colors$Group) <- as.character(c(comp$Test[i], comp$Control[i]))
		if('Fold_Change' %in% names(dif)){
			updown <- subset(dif, p_value <= comp$p_value[i] & q_value <= comp$q_value[i] & abs(log2(Fold_Change)) >= log2(comp$fc[i]))
			oupdown <- updown[order(updown$Fold_Change), ]
			top20 <- rbind(subset(head(oupdown, 20), Fold_Change <= 1/comp$fc[i]), subset(tail(oupdown, 20), Fold_Change >= comp$fc[i]))
		}else if('log2FC' %in% names(dif)){
			updown <- subset(dif, p_value <= comp$p_value[i] & q_value <= comp$q_value[i] & abs(log2FC) >= log2(comp$fc[i]))
			oupdown <- updown[order(updown$log2FC), ]
			top20 <- rbind(subset(head(oupdown, 20), log2FC <= -log2(comp$fc[i])), subset(tail(oupdown, 20), log2FC >= log2(comp$fc[i])))
		}
		heatdf <- updown[,which(colnames(updown) %in% g[,1])]
		heatdf <- as.matrix(heatdf)
		heatdf[is.na(heatdf)] <- 0
		if((max(heatdf)- min(heatdf)) > 30){
			heatdf <- log2(heatdf+1)
		}
		if(! is.null(opt$scale)){
			sc<-opt$scale
		}else if(length(which(anno_col[,1]==vs[1]))>1 && length(which(anno_col[,1]==vs[2]))>1){
			sc<-"row"
		}else if(opt$colum){
			sc<-"column"
		}else{
			sc<-"none"
		}
		if(any(duplicated(g[,1]))){
			anno_col<-NA
		}else{
			rownames(anno_col)<-g[,1]
		}
		if(nrow(heatdf)>1){
			suga <- ga[ga$Group %in% c(comp$Test[i], comp$Control[i]), ]
			oga <- suga[order(suga$Group),]
			heatdf <- heatdf[, as.character(oga$Sample)[oga$Sample %in% colnames(heatdf)]]
			p <- pheatmap(heatdf,show_rownames=show_rowname,annotation_col=anno_col,scale=sc,annotation_legend=T,color = colorpanel(128, colors[1], colors[2], colors[3]),border_color = NA,cluster_cols=opt$clcol,cluster_rows=opt$clrow,fontsize_row=opt$size, annotation_names_col=F, show_colnames=opt$cname, fontsize_col=10, angle_col=opt$angle, silent=T, annotation_colors=anno_colors, main=opt$title)
			ggsave(paste(opt$out,"/heatmap_",vs[1],"_vs_",vs[2],".png",sep=""), plot=as.ggplot(p), width=opt$width, height=opt$height, units='in',dpi=600, bg='white')
			ggsave(paste(opt$out,"/heatmap_",vs[1],"_vs_",vs[2],".pdf",sep=""), plot=as.ggplot(p), width=opt$width, height=opt$height, units='in', bg='white')
			top20 <- heatdf[rownames(top20), as.character(oga$Sample)[oga$Sample %in% colnames(heatdf)]]
			p <- pheatmap(top20,show_rownames=T,annotation_col=anno_col,scale=sc,annotation_legend=T,color = colorpanel(128, colors[1], colors[2], colors[3]),border_color = NA,cluster_cols=opt$clcol,cluster_rows=opt$clrow,silent=T,fontsize_row=8, annotation_names_col=F, show_colnames=opt$cname, fontsize_col=10, angle_col=opt$angle, annotation_colors=anno_colors, , main=opt$title)
			ggsave(paste(opt$out,"/top20_heatmap_",vs[1],"_vs_",vs[2],".png",sep=""), plot=as.ggplot(p), width=opt$width+1, height=opt$height, units='in',dpi=600, bg='white')
			ggsave(paste(opt$out,"/top20_heatmap_",vs[1],"_vs_",vs[2],".pdf",sep=""), plot=as.ggplot(p), width=opt$width+1, height=opt$height, units='in', bg='white')
		}else{
			cat('no enough data',file=paste(opt$out,"/heatmap_",vs[1],"_vs_",vs[2],".log",sep=""))
		}
	}
}else if(file.exists(opt$dir) && is.null(opt$comp)){
	f<-read.table(opt$dir,header=T,sep="\t",check.names=F,row.names=opt$name)
	if(!is.null(opt$pheno)){
		if(length(which((names(f) %in% ga$Sample)))<2){
			f <- data.frame(t(f), check.names=F)
			if(length(which((names(f) %in% ga$Sample)))<2){
				stop("The input file no header or the sample not in group file.")
			}
		}
		ga[,2] <- factor(ga[,2], levels=unique(ga[,2]))
		mfs<-f[,names(f) %in% ga$Sample]
		#mm<-apply(mfs,1,sum)
		#mfs<-mfs[mm>0,]
		if(max(as.matrix(mfs))- min(as.matrix(mfs)) > 30){
			mfs <- log2(mfs+1)
		}
		if(ncol(ga)<3){
			anno_col<-data.frame(ga[,2])
			names(anno_col) <- names(ga)[2]
	
		}else{
			anno_col <- ga[,-1]
		}
		anno_colors <- list(Group=get(paste('pal_', opt$palette, sep=''))()(length(unique(ga[,2])))[as.numeric(unique(ga[,2]))])
		names(anno_colors)[1] <- names(ga)[2]
		names(anno_colors[[names(ga)[2]]]) <- unique(as.character(ga[,2]))
		rownames(anno_col)<-ga$Sample
		ga <- ga[ga$Sample %in% names(mfs),]
		ot <- 'order(ga[, 2]'
		if(ncol(ga) > 2){
			for (i in 3:ncol(ga)){
				ot <- paste(ot, ', ', 'ga[,',i,']', sep='')
			}
		}
		ot <- paste(ot , ')', sep='')
		od <- eval(parse(text = ot))
		oga <- ga[od, ]
		mfs <- mfs[, oga$Sample]
		if(ncol(anno_col)>1){
			anno_col <- anno_col[,rev(colnames(anno_col))]
		}
		if(ncol(anno_col) <2){
			an_col <- F
		}else{
			an_col <- T
		}
	}else{
		mfs <- f
		anno_col <- NA
		anno_colors <- NA

	}
	if(! is.null(opt$scale)){
		sc<-opt$scale
	}else if(!is.na(anno_col) && all(table(anno_col[, ncol(anno_col)])>1)){
		sc<-"row"
	}else if(opt$colum){
		sc<-"column"
	}else{
		sc<-'none'
	}
	if(length(colors)==2){
		clp <-  colorpanel(128, low=colors[1], high=colors[2])
	}else if(length(colors)==3){
		clp <- colorpanel(128, low=colors[1], mid=colors[2], high=colors[3])
	}
	p <- pheatmap(mfs,show_rownames=show_rowname,annotation_col=anno_col,scale=sc,color=clp, border_color = NA,cluster_cols=opt$clcol,cluster_rows=opt$clrow,fontsize_row=opt$size,annotation_legend=T, show_colnames=opt$cname, angle_col=opt$angle, silent=T, annotation_colors=anno_colors, fontsize_col=10, fontsize=10, main=opt$title, legend=opt$legend, annotation_names_col=an_col, clustering_method='average')
	ggsave(paste(opt$out,"heatmap.png",sep=""), plot=as.ggplot(p), width=opt$width, height=opt$height, units='in', dpi=600, bg='white')
	ggsave(paste(opt$out,"heatmap.pdf",sep=""), plot=as.ggplot(p), width=opt$width, height=opt$height, units='in', bg='white')
}
if(file.exists('Rplots.pdf')){
	file.remove("Rplots.pdf")
}
