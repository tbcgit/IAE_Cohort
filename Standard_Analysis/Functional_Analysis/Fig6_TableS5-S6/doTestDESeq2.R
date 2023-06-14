#!/usr/bin/env Rscript
###############

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
myfactors <- as.character(parameters[3])
mylevels <- as.character(parameters[4])
study <- as.character(parameters[5])
doFilt <- as.logical(as.character(parameters[6]))
mypct <- as.numeric(as.character(parameters[7]))
myfac <- as.numeric(as.character(parameters[8]))
descr.column <- as.character(parameters[9])
mypaired <- as.logical(as.character(parameters[10]))
mynorm <- as.logical(as.character(parameters[11]))
myfitType <- as.character(parameters[12])

source("myfunctions.downstream.R")
library(DESeq2)

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
DESCRIPTION <- dat[, descr.column]; names(DESCRIPTION) <- rownames(dat) <- dat[, 1]

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]
mysamples <- intersect(rownames(meta), colnames(dat))
dat <- as.matrix(dat[, mysamples])
meta <- meta[mysamples, ]

myfactors.l <- strsplit(myfactors, split=",")
myfactor1 <- myfactors.l[[1]][1]
if(length(myfactors.l[[1]])==2) myfactor2 <- myfactors.l[[1]][2] else myfactor2 <- ""

if(myfactor2!="") GROUPS <- paste(meta[, myfactor1], meta[, myfactor2], sep=".") else GROUPS <- as.character(meta[, myfactor1])
names(GROUPS) <- mysamples

mylevels <- unlist(strsplit(mylevels, split=","))

mylevels <- mylevels[which(is.element(mylevels, unique(GROUPS)))]

ind <- which(is.element(GROUPS, mylevels))

dat <- dat[, ind]
GROUPS <- GROUPS[ind]
meta <- meta[ind, ]

# Esto deja los datos exactamente igual si la matriz era una matriz de counts relativos (en porcentajes, escala 0-100, 0-1, etc).
# Eliminamos genes que tienen un nivel de expresión asignado demasiado alto.
# Dejamos a cero los que tienen expresión demasiado baja.
# Por lo visto DESeq2 da problemas con números muy grandes y/o con una distribución de datos aberrantemente polarizada.
mymin <- min(dat[which(dat[, 1]!=0), 1])
mymax <- max(dat[, 1])
mysum <- sum(dat[, 1])
lo <- 1000
if(mymin!=1&mysum!=100&mysum!=1&(mymax/mymin)>100000){
	aberrant <- c()
	#mymins <- apply(dat, MARGIN=2, FUN=function(x) min(x[which(x!=0)]))
	#mycutdown <- max(mymins)*5
	#dat[which(dat<mycutdown)] <- 0
	for(i in 1:ncol(dat)){
		mycutdown <- min(dat[which(dat[, i]!=0), i])*5
		dat[which(dat[, i]<mycutdown), i] <- 0
		if(nrow(dat)>lo) qs <- quantile(dat[, i], probs=seq(from=0, to=1, length.out=lo)) else qs <- quantile(dat[, i], probs=seq(from=0, to=1, length.out=nrow(dat)-10))
		mycutup <- qs[length(qs)-1]
		# seleccionamos como aberrantes el 0.1% más grandes
		aberrant <- c(aberrant, rownames(dat)[which(dat[, i]>mycutup)])
	}
	aberrant <- unique(aberrant)
	dat <- dat[setdiff(rownames(dat), aberrant), ]
	#for(i in 1:ncol(dat)) dat[, i] <- round(dat[, i]/min(dat[which(dat[, i]!=0), i]))
	write(aberrant, file="aberrant.features.txt", sep="\n")
}

tem <- apply(dat, MARGIN=2, sum)
tem <- abs(tem-tem[1])
if(all(tem<0.0001)){
	# Cuando todas las muestras suman lo mismo (o muy parecido) se supone que en la situación previa
	# de freqs absolutas (antes de pasar a relativas) cada muestra tenía diferente número de lecturas.
	# Esta transformación es correcta si en todas las muestras había algún 1 en los datos originales en freqs absolutas.
	for(i in 1:ncol(dat)) dat[, i] <- round(dat[, i]/min(dat[which(dat[, i]!=0), i]))
} else {
	# Si la entrada es una matriz de frecuencias absolutas con algún 1, esto no cambia los datos.
	# En caso de que no todas las muestras sumen lo mismo (o muy parecido) es mejor dejar los datos igual pero en escala absoluta.
	dat <- round(dat/min(dat[which(dat!=0)]))
}

cond <- apply(dat, MARGIN=1, sum)
dat <- dat[which(cond>0), ]

if(any(as.vector(as.matrix(dat))<0)) doCCA <- FALSE else doCCA <- TRUE

if(mynorm==TRUE) datnorm <- prop.table(as.matrix(dat), margin=2)*100 else datnorm <- dat

if(doFilt==TRUE){
	mycut <- min(abs(datnorm[which(datnorm!=0)]))*myfac
	good <- doFiltering(datnorm, cutoff=mycut, group=GROUPS, mypct=mypct)
	print(length(good))
	dat <- round(dat[good, ], digits=5)
	datnorm <- round(datnorm[good, ], digits=5)
	mydata <- as.data.frame(dat)
	mydata$ID <- rownames(mydata)
	mydata <- mydata[, c("ID", colnames(dat))]
	write.table(mydata, file=paste("filtered.freqs", counts.file, sep="."), sep="\t", row.names=FALSE, quote=FALSE)
	cond <- apply(dat, MARGIN=2, sum)
	dat <- dat[, which(cond>0)]; datnorm <- datnorm[, which(cond>0)]
	meta <- meta[which(cond>0), ]; GROUPS <- GROUPS[which(cond>0)]
}

if(length(mylevels)>2){
	if(doCCA==TRUE) correspondence.analysis(MATRIX=datnorm,
					GROUPS=GROUPS, FILE=paste("adonis.CCA", "all.groups", study, "pdf", sep="."), plot.sites.labels=FALSE)
}

mycombn <- combn(mylevels,2)

for(i in 1:ncol(mycombn)){

	g1 <- mycombn[1,i]
    	g2 <- mycombn[2,i]
	i1 <- which(GROUPS==g1)
	i2 <- which(GROUPS==g2)

	if(mypaired==TRUE){
		names(i1) <- meta$donor[i1]
		names(i2) <- meta$donor[i2]
		mydonors <- intersect(names(i1), names(i2))
		i1 <- i1[mydonors]; i2 <- i2[mydonors]
	}

	ind <- c(i1, i2)

	tag <- paste(g1, "vs", g2, sep=".")
	
	if(doCCA==TRUE) correspondence.analysis(MATRIX=datnorm[, ind],
					GROUPS=GROUPS[ind], FILE=paste("adonis.CCA", tag, study, "pdf", sep="."), plot.sites.labels=FALSE)
	
	mytb <- differential.tests(MATRIX=dat[, ind], GROUPS=GROUPS[ind], GROUP1=g1, GROUP2=g2, TEST="DESeq2.test",
						PAIRED=mypaired, CUTOFF=0.1, COLNAME="feature", FILE=paste("DESeq2.test", tag, study, sep="."),
						DESCRIPTION=DESCRIPTION, myfitType=myfitType, meta=meta[ind, ], datnorm=datnorm[, ind])

	myfeatures <- as.character(mytb[which(mytb$p.value<=0.05), "feature"])
	tem <- as.data.frame(dat[myfeatures, ind]); tem$ID <- rownames(tem); tem <- tem[, c("ID", setdiff(colnames(tem), "ID"))]
	write.table(tem, file=paste("significant.rel.freqs", tag, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE) 
	
	if(mypaired==TRUE){
		pdf(file=paste("DESeq.test.paired.with.lines", tag, study, "pdf", sep="."))
		for(myfeature in myfeatures){
			plotPaired(dat=datnorm[, ind], feature=myfeature, group=GROUPS[ind],
					group.x=g1, group.y=g2, individual.factor=meta$donor[ind], DESCRIPTION=DESCRIPTION)
		}
		dev.off()
	}

	print(tag)
}




