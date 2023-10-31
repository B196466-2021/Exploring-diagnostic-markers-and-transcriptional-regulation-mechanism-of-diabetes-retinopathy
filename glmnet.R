#! /usr/bin/Rscript

usage<-function(spec){
	cat("This script is used to\n",
getopt(spec,usage=TRUE),
"Options:
	-f, --file	input matrix (sample in each column and variable in each row).
	-r, --resp	response variable file (see details at 'https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html' for each response type).
	-l, --lvf	list file of variable will be used (e.g. gene).
	-i, --image	R image of fited results (only for exporting fited results).
	-a, --ratio	the ratio to split the sample into train and test group (default 0.5, 1 means do not split).
	-d, --rds	rds file of cv.glmnet which was used to validate the fit result (without -l)
	-t, --type	response type: 'gaussian', 'binomial', 'poisson', 'multinomial', 'cox', 'mgaussian' or else a 'glm()' family (default 'binomial').
	-e, --corr	correlation threshold value used to filter genes (defaut 0.75).
	-o, --out	prefix of out files(default ./).
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
	'resp','r',1,'character',
	'lvf','l',1,'character',
	'ratio','a',1,'double',
	'image','i',1,'character',
	'rds','d',1,'character',
	'type','t',1,'character',
	'corr','e',1,'double',
	'out','o',1,'character',
	'version','v',0,'logical',
	'help','h',0,'logical'
	),byrow=TRUE, ncol=4)

opt=getopt(spec)

if (!is.null(opt$version)) {version()}
if (is.null(opt$image) && (is.null(opt$file) || is.null(opt$resp)) || !is.null(opt$help)) {usage(spec)}
if (is.null(opt$out)) {opt$out='./'}
if (is.null(opt$ratio)) {opt$ratio=0.5}
if (is.null(opt$corr)) {opt$corr=0.75}
if (is.null(opt$type)) {opt$type='binomial'}
dir <- gsub("[^/]*$", '', opt$out, perl=T)
if(! dir.exists(dir) && dir != ''){
	dir.create(dir, recursive=T)
}
##########################
#define functions
rocfun <- function(tm, dsr){
	ucox <- coxph(Surv(time, status) ~ risk, data = dsr)
    dsr$lp <- predict(ucox, type = "lp")
	roc <- survivalROC(Stime = dsr$time, status = dsr$status, marker = dsr$lp, predict.time = tm, span = 0.25*nrow(dsr)^(-0.20))
	auct <- data.frame(time=tm, auc=roc$AUC)
	return(auct)
}


##########################
#main programe
library(glmnet)
library(Hmisc)
library(pROC)
library(caret)
library(ggsci)
#import file
if(!is.null(opt$rds) && file.exists(opt$rds)){
	cvfit <- readRDS(opt$rds)
	val <- rownames(cvfit$glmnet.fit$beta)
}else if(!is.null(opt$lvf) && file.exists(opt$lvf)){
	val <- readLines(opt$lvf)
}
pred.type <- list(binomial='class', gaussian='link', mgaussian='link', cox='response', poisson='response')

if(!is.null(opt$file) && exists('val')){
	mda <- read.table(opt$file, sep='\t', header=T, row.names=1, check.names=F, quote="")
	sdat <- t(mda[val,])
	colnames(sdat) <- val
	sdat[is.na(sdat)] <- 0
	if(is.null(opt$resp)){
		stop('no response file with -r or --resp.')
	}
	resp <- read.table(opt$resp, sep='\t', header=T, row.names=1, check.names=F)
	sdats <- sdat[rownames(sdat) %in% rownames(resp), ]
	old.name <- colnames(sdats)
	new.name <- gsub('-', '_', colnames(sdats))
	colnames(sdats) <- new.name
	respo <- matrix(resp[rownames(sdats),])
	rownames(respo) <- rownames(sdats)
	colnames(respo) <- colnames(resp)
}
if(is.null(opt$rds) && is.null(opt$image) && !is.null(opt$file)){
	#filter high correlated features
	descrCor <- cor(sdats)
	highlyCorDescr <- findCorrelation(descrCor, cutoff=opt$corr)
	gene <- colnames(descrCor)[-highlyCorDescr]
	print(paste(length(gene), 'gene remained after filter by correlation with', opt$corr))
	if(length(gene) < 5){
		stop(paste('remained genes is less than 5, increase the correlation threshold.'))
	}
	sdats <- sdats[, gene]
	#fit
	if(opt$ratio ==1){
		fit <- glmnet(sdats, respo, family=opt$type)
		cvfit <- cv.glmnet(sdats, respo, family=opt$type)
		pred.all <- predict(cvfit, newx=sdats, s = "lambda.min", type=pred.type[[opt$type]])
		rrsd <- merge(resp, pred.all, by=0)
		names(rrsd) <- c('Sample', names(resp), 'Predictor')
		roc.overall <- roc(response=respo, predictor=ordered(pred.all), direction='<')
		if(as.numeric(roc.overall$auc) < 0.5){
			roc.overall <- roc(response=respo, predictor=ordered(pred.all), direction='>')
		}
		pd <- data.frame(x=roc.overall$specificities, y=roc.overall$sensitivities, AUC=paste('AUC (', format(as.numeric(roc.overall$auc), digits=3), ')', sep=''))
	}else if(opt$ratio > 0 && opt$ratio <1){
		auc.train <- auc.test <- 0.65
		time1 <- proc.time()[3]
		time2 <- proc.time()[3]
		tn <- 1
		auco <- data.frame()
		tmd <- paste(opt$out, '/tmp/', sep='')
		if(!dir.exists(tmd)){
			dir.create(tmd, recursive=T)
		}
		writeLines('train_id\tauc.train\tauc.test\tauc.all', paste(tmd, '/auc.txt', sep=''))
		while((as.numeric(time2 - time1) < 230400 || auc.train < 0.7 || auc.test < 0.7 || auc.all < 0.7) && as.numeric(time2 - time1) < 430000){
			if(auc.train==1 && auc.test==1){
				break
			}
			if(opt$type=='binomial'){
				trainIndex <- createDataPartition(factor(respo[,1]), p=opt$ratio, list=FALSE, times=1)
				trs <- rownames(respo)[trainIndex]
			}else{
				trainIndex <- sample(rownames(respo), size=round(nrow(respo)*opt$ratio))
			}
			trd <- sdats[trs, ]
			ted <- sdats[!rownames(sdats) %in% trs, ]
			train.res <- respo[rownames(trd),]
			test.res <- respo[rownames(ted),]
			cvfit <- cv.glmnet(x=trd, y=train.res, family=opt$type)
			risk.train <- predict(cvfit, newx=trd, s = "lambda.min", type=pred.type[[opt$type]])
			risk.test <- predict(cvfit, newx=ted, s = "lambda.min", type=pred.type[[opt$type]])
			risk.all <- predict(cvfit, newx=sdats, s = "lambda.min", type=pred.type[[opt$type]])
			roc.tr <- roc(response=train.res, predictor=ordered(risk.train))
			auc.tr <- as.numeric(roc.tr$auc)
			roc.te <- roc(response=test.res, predictor=ordered(risk.test))
			auc.te <- as.numeric(roc.te$auc)
			roc.all <- roc(response=respo, predictor=ordered(risk.all))
			auc.all <- as.numeric(roc.all$auc)
			#if(auc.tr < 0.5 && auc.te < 0.5){
			#	roc.tr <- roc(response=train.res, predictor=ordered(risk.train), direction='>')
			#	roc.te <- roc(response=test.res, predictor=ordered(risk.test), direction='>')
			#	roc.all <- roc(response=respo, predictor=ordered(risk.all), direction='>')
			#	auc.tr <- as.numeric(roc.tr$auc)
			#	auc.te <- as.numeric(roc.te$auc)
			#	auc.all <- as.numeric(roc.all$auc)
			#}
			if(auc.tr + auc.te > auc.train + auc.test && abs(auc.tr - auc.te) < 0.05){
				auc.train <- auc.tr
				auc.test <- auc.te
				roc.train <- roc.tr
				roc.test <- roc.te
				roc.overall <- roc.all
				pred.train <- risk.train
				pred.test <- risk.test
				pred.all <- risk.all
				auc.sub <- data.frame(train_id=tn, auc.train=auc.train, auc.test=auc.test, auc.all=auc.all)
				write.table(auc.sub, paste(tmd, '/auc.txt', sep=''), sep='\t', row.names=F, quote=F, append=T, col.names=F)
				auco <- rbind(auco, auc.sub)
				cvfit.used <- cvfit
				fit <- glmnet(trd, train.res, family=opt$type)
				save.image(file = paste(tmd, '/train_', tn, '.RData', sep=''))
				tn <- tn + 1
			}
			time2 <- proc.time()[3]
		}
	}
}else if(!is.null(opt$image)){
	load(opt$image)
	if(!exists('fit')){
		stop('no train result is available.')
	}
}
if(is.null(opt$rds) && exists('cvfit.used') && opt$ratio !=1){
	cvfit <- cvfit.used
	rrsd <- merge(resp, pred.train, by=0)
	names(rrsd) <- c('Sample', names(resp), 'Predictor')
	write.table(rrsd, file=paste(opt$out, 'glmnet_train_predictor.txt', sep=''), quote=F, sep='\t', row.names=F)
	rrsd <- merge(resp, pred.test, by=0)
	names(rrsd) <- c('Sample', names(resp), 'Predictor')
	write.table(rrsd, file=paste(opt$out, 'glmnet_test_predictor.txt', sep=''), quote=F, sep='\t', row.names=F)
	pd <- data.frame(x=c(roc.train$specificities, roc.test$specificities, roc.overall$specificities), 
					 y=c(roc.train$sensitivities, roc.test$sensitivities, roc.overall$sensitivities), 
					 AUC=c(rep(paste('Train AUC (', format(as.numeric(roc.train$auc), digits=3), ')', sep=''), times=length(roc.train$specificities)),
						   rep(paste('Test AUC (', format(as.numeric(roc.test$auc), digits=3), ')', sep=''), times=length(roc.test$specificities)),
						   rep(paste('All AUC (', format(as.numeric(roc.overall$auc), digits=3), ')', sep=''), times=length(roc.overall$specificities))))
	pd$AUC <- factor(pd$AUC, levels=sort(unique(pd$AUC), decreasing=T))
}

if(is.null(opt$rds) && exists('pd')){
	pd <- pd[order(1-pd$x, pd$y),]
	p <- ggplot(pd, aes(x=1-x, y=y, color=AUC)) + geom_line(size=1)
	p <- p + geom_segment(x=0, xend=1, y=0, yend=1, size=0.5, linetype=2, color='grey40')
	p <- p + labs(title=opt$title, x='1-Specificity', y='Sensitivity') + scale_color_npg() + theme_bw()
	p <- p + theme(axis.title=element_text(size=20), axis.text=element_text(size=15, color='black'), legend.text=element_text(size=12), 
				   legend.position=c(0.65, 0.2), plot.title=element_text(size=15, hjust=0.5), legend.title=element_blank())
	ggsave(paste(opt$out, 'roc.png', sep=''), plot=p, width=4, height=4, units='in', dpi=600)
	ggsave(paste(opt$out, 'roc.pdf', sep=''), plot=p, width=4, height=4, units='in')
	rrsd <- merge(resp, pred.all, by=0)
	names(rrsd) <- c('Sample', names(resp), 'Predictor')
	if(pred.type[[opt$type]]=='class'){
		pb <- predict(cvfit, newx=sdats, s = "lambda.min", type='response')
		rrsd <-  merge(rrsd, pb, by.x=1, by.y=0)
		names(rrsd)[ncol(rrsd)] <- 'Probability'
	}
	write.table(rrsd, file=paste(opt$out, 'glmnet_all_predictor.txt', sep=''), quote=F, sep='\t', row.names=F)
	pdf(paste(opt$out, 'lasso_fit.pdf', sep=''), width=5, height=5)
	plot(fit, xvar="lambda", label=TRUE)
	abline(v=log(cvfit$lambda.min), col='red')
	text(x=log(cvfit$lambda.min), y=min(fit$beta), labels=paste('lambda = ', round(cvfit$lambda.min, 4), sep=''), col='red')
	dev.off()

	png(paste(opt$out, 'lasso_fit.png', sep=''), res=600, units='in', width=5, height=5)
	plot(fit, xvar="lambda", label=TRUE)
	abline(v=log(cvfit$lambda.min), col='red')
	text(x=log(cvfit$lambda.min), y=min(fit$beta), labels=paste('lambda = ', round(cvfit$lambda.min, 4), sep=''), col='red')
	dev.off()

	png(paste(opt$out, 'lasso_cross_validation.png', sep=''), res=600, units='in', width=5, height=5)
	plot(cvfit)
	dev.off()
		pdf(paste(opt$out, 'lasso_cross_validation.pdf', sep=''), width=5, height=5)
	plot(cvfit)
	dev.off()
		coef.min = as.matrix(coef(cvfit, s = "lambda.min"))
	active.min = which(coef.min != 0)
	od <- data.frame(Gene=row.names(coef.min)[active.min], Coef=coef.min[active.min])
	if(nrow(od) < 1 ){
		stop('No gene retained.')
	}
	lco <- data.frame(variable=row.names(coef.min), Coef=coef.min[,1])
	lco$variable <- translate(lco$variable, new.name, old.name)
	write.table(lco, file=paste(opt$out, 'glmnet_coefficient.txt', sep=''), quote=F, sep='\t', row.names=F)
	rownames(cvfit$glmnet.fit$beta) <- translate(rownames(cvfit$glmnet.fit$beta), new.name, old.name)
	saveRDS(cvfit, file=paste(opt$out, 'cvfit_result.rds', sep=''))
}else if(!is.null(opt$rds) && exists('cvfit')){
	#validate
	pred.all <- predict(cvfit, newx=sdats, s = "lambda.min", type=pred.type[[opt$type]])
	rrsd <- merge(resp, pred.all, by=0)
	names(rrsd) <- c('Sample', names(resp), 'Predictor')
	if(pred.type[[opt$type]]=='class'){
		pb <- predict(cvfit, newx=sdats, s = "lambda.min", type='response')
		rrsd <-  merge(rrsd, pb, by.x=1, by.y=0)
		names(rrsd)[ncol(rrsd)] <- 'Probability'
	}
	write.table(rrsd, file=paste(opt$out, 'glmnet_all_predictor.txt', sep=''), quote=F, sep='\t', row.names=F)
	roc.overall <- roc(response=respo, predictor=ordered(pred.all))
	#if(as.numeric(roc.overall$auc) < 0.5){
	#	roc.overall <- roc(response=respo, predictor=ordered(pred.all), direction='>')
	#}
	pd <- data.frame(x=roc.overall$specificities, y=roc.overall$sensitivities, AUC=paste('AUC (', format(as.numeric(roc.overall$auc), digits=3), ')', sep=''))
	pd <- pd[order(1-pd$x, pd$y),]
	p <- ggplot(pd, aes(x=1-x, y=y, color=AUC)) + geom_line(size=1)
	p <- p + geom_segment(x=0, xend=1, y=0, yend=1, size=0.5, linetype=2, color='grey40')
	p <- p + labs(title=opt$title, x='1-Specificity', y='Sensitivity') + scale_color_npg() + theme_bw()
	p <- p + theme(axis.title=element_text(size=20), axis.text=element_text(size=15, color='black'), legend.text=element_text(size=12),
				   legend.position=c(0.65, 0.2), plot.title=element_text(size=15, hjust=0.5), legend.title=element_blank())
	ggsave(paste(opt$out, 'roc.png', sep=''), plot=p, width=4, height=4, units='in', dpi=600)
	ggsave(paste(opt$out, 'roc.pdf', sep=''), plot=p, width=4, height=4, units='in')
	pr <- merge(pred.all, respo, by=0)
	names(pr)[1:2] <- c('Sample', 'Probabilities')
	write.table(pr, file=paste(opt$out, 'cvfit_validate_predict_result.txt', sep=''), sep='\t', quote=F, row.names=F)
}


