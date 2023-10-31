#! /usr/bin/Rscript

##########################
#get options

library(getopt)
opt = getopt(matrix(c(
'dir','d',1,'character',
'comp','c',1,'character',
'type','t',1,'character',
'color','r',1,'character',
'lid','i',1,'character',
'pval','l',1,'double',
'qval','q',1,'double',
'fold','f',1,'double',
'xlab','x',1,'character',
'ylim','y',1,'double',
'out','o',1,'character',
'version','v',0,'logical',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));

usage<-function(){
	cat("This script is used to plot differential expression.
Usage   Rscript",get_Rscript_filename(),"[options]
Options:
	-d, --dir	directory congtaining all differential result file of each compare, or differential result file of one compare.
	-c, --comp	file contain all compare informations.
	-r, --color	color for up and down regulating(default brown3,blue3).
	-i, --lid	list file of id or the number of column will be displayed as text.
	-l, --pval	threshold value of p_value (defult 0.05).
	-q, --qval	threshold value of q_value.
	-f, --fold	the fold of up and down regulate(defult 1.5).
	-x, --xlab	the label of x for vocano plot (defult log2(Fold Change)).
	-o, --out	output directory(default the program runing directory).
	-y, --ylim	limit of y axis for volcano plot (default max vlaue).
	-v, --version	display version informations.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1);
}

if (!is.null(opt$version)) {stop(version())}
if (is.null(opt$dir) || !is.null(opt$help)){
	stop(usage())
}

if (is.null(opt$type)) {opt$type="protein_coding"}
if (is.null(opt$col)) {opt$col="brown3,blue3"}
if (is.null(opt$out)) {opt$out="."}
if (is.null(opt$pval)) {opt$pval=0.05}
if (is.null(opt$fold)) {opt$fold=1.5}
if (is.null(opt$xlab)) {opt$xlab='log2(Fold Change)'}
##########################
#subfunctions
scatter<-function(x,test=1,control=2,times=1.5,xlab="control",ylab="test",labe=NULL,log2FC=4){
	if(max(data.frame(x[,control],x[,test]))>25){
		x[,control]<-log2(x[,control]+1)
		x[,test]<-log2(x[,test]+1)
		diagonal<-FALSE
	}else{
		diagonal<-TRUE
	}
	ud<-x[,log2FC]
	ud[which(ud>=log2(times))]<-100
	ud[which(ud<=log2(1/times))]<-"down"
	ud[which(ud=="100")]<-"up"
	ud[grep("up|down",ud,invert=T)]<-"nde"
	sca<-data.frame(test=x[,test],Control=x[,control],ud,x[,1])
	corr<-cor(sca[,1:2])
	up<-paste("up (",length(grep("up",ud)),")",sep="")
	down<-paste("down (",length(grep("down",ud)),")",sep="")
	nde<-paste("middle (",length(grep("nde",ud)),")",sep="")
	p<-ggplot(sca,aes(x=Control,y=test,colour=ud))
	p<-p+geom_point(shape=20,size=3)
	p<-p+scale_colour_manual(values=c("up"=color[1],"nde"="grey60","down"=color[2]),breaks = c("up","nde","down"),labels=c(up,nde,down))
	p<-p+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
	if(diagonal){
		p<-p+geom_abline(intercept = c(log2(times),-log2(times)), slope = 1,colour="grey45",linetype=2,size=0.3)
	}
	p<-p+labs(x=xlab,y=ylab,title=paste(ylab,"vs",xlab,by=" "))
	p<-p+guides(col = guide_legend(title=paste("pearson correlation:",round(corr[1,2],3),sep=" "),title.hjust=0.3,title.theme = element_text(size=12,angle = 0), ncol =1, override.aes=list(size=5)))
	p<-p+theme(legend.position=c(0.25,0.85),legend.text = element_text(size = 15))+ylim(min(c(sca$test, sca$Control)), round(max(c(sca$test, sca$Control))*1.2))+xlim(min(c(sca$test, sca$Control)), round(max(c(sca$test, sca$Control))*1.2))
	p<-p+theme(plot.margin=margin(15, 40, 10, 10), axis.title=element_text(size=20), axis.text=element_text(size=15, colour='black'), plot.title=element_text(size=20))
	if(!is.null(labe)){
		lab<-merge(sca,labe, by.x=4, by.y=1)
		lab$ud<-gsub("up",color[1],lab$ud)
		lab$ud<-gsub("down",color[2],lab$ud)
		lab$ud<-gsub("nde","grey60",lab$ud)
		p<-p+geom_text_repel(data=lab,aes(x=Control,y=test,label=lab[,5]),colour="black",inherit.aes = FALSE,show.legend=F,size=3, min.segment.length=0)
	}
	return(p)
}

volcano<-function(x,log2FC=1,value=2,level=0.05,times=1.5,title="test vs control",labe=NULL){
	du2<-x[,log2FC]
	du2[which(x[,value]>level)]<-0
	du2[which(du2>=log2(times))]<-100
	du2[which(du2<=(-log2(times)))]<-"down"
	du2[which(du2=="100")]<-"up"
	du2[which(!du2 %in% c("up","down"))]<-"nde"
	vol2<-data.frame(x[,log2FC],-log10(x[,value]),du2,id=x[,1])
	upl<-paste("up (",length(grep("up",du2)),")",sep="")
	downl<-paste("down (",length(grep("down",du2)),")",sep="")
	ndel<-paste("middle (",length(grep("nde",du2)),")",sep="")
	p<-ggplot(vol2,aes(x=vol2[,1],y=vol2[,2],colour=du2))
	p<-p+geom_point(shape=20,size=3)
	p<-p+scale_colour_manual(values=c("up"=color[1],"nde"="grey60","down"=color[2]),breaks = c("up","nde","down"),labels=c(upl,ndel,downl))
	p<-p+geom_hline(yintercept = -log10(level),colour="green4",linetype=2,size=0.3)
	p<-p+geom_vline(xintercept = c(-log2(times),log2(times)),colour="green4",linetype=2,size=0.3)
	p<-p+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
	p<-p+labs(x=opt$xlab,y=paste("-log10(",names(x)[value],")",sep=""),title=title)
	p<-p+theme(legend.position="top",legend.text = element_text(size = 15))
	p<-p+guides(col = guide_legend(nrow = 1,title=NULL, override.aes=list(size=5)))
	p<-p+theme(plot.margin=margin(15, 40, 10, 10), axis.title=element_text(size=20), axis.text=element_text(size=15, colour='black'), plot.title=element_text(size=20))
	p<-p+xlim(-(max(abs(vol2[,1]))),max(abs(vol2[,1])))
	if(!is.null(labe)){
		lab<-merge(vol2,labe, by.x=4, by.y=1)
		lab$du2<-gsub("up",color[1],lab$du2)
		lab$du2<-gsub("down",color[2],lab$du2)
		lab$du2<-gsub("nde","grey60",lab$du2)
		p<-p+geom_text_repel(data=lab,aes(x=lab[,2],y=lab[,3],label=lab[,5]),colour="black",inherit.aes = FALSE,show.legend=F,size=3, min.segment.length=0)
	}
	return(p)
}

#######################
#processing data

library(ggplot2)
library(ggrepel)
library(dplyr)
color<-strsplit(opt$col, ',')[[1]]
if(is.null(opt$lid)){
	labe<-NULL
}else if(file.exists(opt$lid)){
	labe<-read.table(opt$lid,header=F,sep="\t")
	if(ncol(labe)==1){
		labe<-data.frame(labe,labe)
	}
}

lncRNA<-c("processed_transcript","lincRNA","3prime_overlapping_ncrna","antisense","non_coding","sense_intronic","sense_overlapping","TEC","known_ncrna","macro_lncRNA","bidirectional_promoter_lncrna")
sncRNA<-c("snRNA","snoRNA","rRNA","Mt_tRNA","Mt_rRNA","misc_RNA","miRNA","ribozyme","sRNA","scaRNA","vaultRNA")
otp<-c("Gene_Type","Trans_Type")

if(!is.null(opt$comp)){
	comp<-read.table(opt$comp,header=T,sep="\t",colClasses = c(rep("character", 2), rep("numeric", 3)))
	for(i in 1:nrow(comp)){
		if(dir.exists(opt$dir)){
			filep<-paste(opt$dir, '/', comp$Test[i], '_vs_', comp$Control[i], '.txt', sep='')
		}else{
			filep<-opt$dir
		}
		dif<-read.table(filep, header=T, sep='\t', check.names=F, stringsAsFactors=F, quote="")
		tm<-otp[otp %in% names(dif)]
		if(length(tm)==0){
			tm<-"Gene_Type"
		}
		if(tm %in% names(dif)){
			if(opt$type=="all"){
					type<-unique(dif[,eval(tm)])
			}else{
				type<-strsplit(opt$type,",")[[1]]
				if("lncRNA" %in% type){
					type<-c(type,lncRNA)
				}
				if("sncRNA" %in% type){
					type<-c(type,sncRNA)
				}
			}
		}
		if(tm %in% names(dif)){
			m3<-dif[m[,eval(tm)] %in% type,]
			if(nrow(m2)==0){
				stop(paste("no Gene_Type name(s): ",paste(type,collapse=" "),sep=""))
			}
			nit<-type[which(!type %in% unique(m2[,eval(tm)]))]
			if(length(nit) != 0){
				cat(paste("Warning: no Gene_Type name(s): ",paste(nit,collapse=" "),"\n",sep=""))
			}
		}else{
			m3<-dif
		}
		if(!is.null(opt$lid) && !file.exists(opt$lid)){
			uplab <- m3 %>% filter(q_value < comp$q_value[i] & p_value < comp$p_value[i])  %>% top_n(n=5, wt=log2FC)
			downlab <- m3 %>% filter(q_value < comp$q_value[i] & p_value < comp$p_value[i])  %>% top_n(n=5, wt=-log2FC)
			labe <- rbind(uplab, downlab)[, c(1, as.numeric(opt$lid))]
			names(labe) <- c('id', 'label')
		}
		#if(all(c(comp$Test[i], comp$Control[i]) %in% names(m3))){
		#	p<-scatter(m3,test=which(names(m3)==comp$Test[i]),control=which(names(m3)==comp$Control[i]),times=comp[i,3],xlab=comp[i,2],ylab=comp[i,1],labe=labe,log2FC=which(names(m3)=="log2FC"))
		#	ggsave(paste(opt$out,"/scatter_",comp[i,1],"_vs_",comp[i,2],".pdf",sep=""),plot=p,device="pdf",width=6,height=6,units="in")
		#	ggsave(paste(opt$out,"/scatter_",comp[i,1],"_vs_",comp[i,2],".png",sep=""),plot=p,device="png",width=6,height=6,units="in",dpi=600)
		#}else{
		#	cat("Warning: no FPKM value in Test_value or Control_value\n")
		#}

		#if(all(m3$p_value==1) && all(m3$q_value==1)){
		if(comp$p_value[i]==1 && comp$q_value[i]==1){
			next
		}

		if(((comp$q_value[i]==1 || comp$p_value[i] < comp$q_value[i] || !"q_value" %in% names(comp)) && all(c("log2FC","p_value") %in% names(m3))) || !"q_value" %in% names(m3)){
			value<-"p_value"
			level<-comp$p_value[i]
		}else if(((comp$p_value[i]==1 || comp$q_value[i] <= comp$p_value[i] || !"p_value" %in% names(comp)) && all(c("log2FC","q_value") %in% names(m3))) || !"p_value" %in% names(m3)){
			value<-"q_value"
			level<-comp$q_value[i]
		}else{
			next
		}
		if(! is.null(opt$ylim)){
			m3[, value][-log10(m3[, value]) > opt$ylim] <- 10**(-opt$ylim)
		}
		p<-volcano(m3,log2FC=which(names(m3)=="log2FC"),value=which(names(m3)==value),level=level,times=comp[i,3],title=paste(comp[i,1],"vs",comp[i,2]),labe=labe)
		ggsave(paste(opt$out,"/volcano_",comp[i,1],"_vs_",comp[i,2],".pdf",sep=""),plot=p,device="pdf",width=6,height=6.5,units="in")
		ggsave(paste(opt$out,"/volcano_",comp[i,1],"_vs_",comp[i,2],".png",sep=""),plot=p,device="png",width=6,height=6.5,units="in",dpi=600)
	}
}else{
	valnam<-names(dif)[max(which(names(dif) %in% c('log2FC',  'p_value', 'q_value)')))+c(1, 2)]
	if(length(valnam) != 2){
		stop(paste("The number of element isn't two which contain '_FPKM' in the header of ",opt$file,"\n",sep=""))
	}
	lab<-valnam
	if(length(lab)==2){
		p<-scatter(dif,test=which(names(dif)==valnam[1]),control=which(names(dif)==valnam[2]),times=opt$fold,xlab=lab[2],ylab=lab[1],labe=labe)
		ggsave(paste(opt$out,"/scatter_",lab[1],"_vs_",lab[2],".pdf",sep=""),plot=p,device="pdf",width=6,height=6,units="in")
		ggsave(paste(opt$out,"/scatter_",lab[1],"_vs_",lab[2],".png",sep=""),plot=p,device="png",width=6,height=6,units="in",dpi=600)
	}else{
		stop(paste("The format of ",opt$file," isn't right\n",sep=""))
	}

	if(all(c("log2FC","p_value", "q_value") %in% names(dif))){
		if((!is.null(opt$qval))){
			vn<-'q_value'
		}else{
			vn<-'p_value'
		}
	}else if(all(c("log2FC","p_value") %in% names(dif))){
		vn<-'p_value'
	}
	if(exists('vn')){
		p<-volcano(dif,log2FC=which(names(dif)=="log2FC"),value=which(names(dif)==vn),level=opt$value,times=opt$fold,title=paste(lab[1],"vs",lab[2],sep=" "),labe=labe)
		ggsave(paste(opt$out,"/volcano_",lab[1],"_vs_",lab[2],".pdf",sep=""),plot=p,device="pdf",width=6,height=6.5,units="in")
		ggsave(paste(opt$out,"/volcano_",lab[1],"_vs_",lab[2],".png",sep=""),plot=p,device="png",width=6,height=6.5,units="in",dpi=600)
	}else{
		cat(paste("Warning: log2FC and p_value or q_value aren'n all in the header of ",opt$file,"\n",sep=""))
	}
}
if(file.exists('Rplots.pdf')){
	file.remove('Rplots.pdf')
}
