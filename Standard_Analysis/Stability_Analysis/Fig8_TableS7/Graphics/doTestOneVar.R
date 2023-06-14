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

mynorm <- FALSE

source("myfunctions.downstream.R")

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

if(mynorm==TRUE) dat <- prop.table(as.matrix(dat), margin=2)*100

dat <- t(as.matrix(dat))
rownames(dat) <- names(DESCRIPTION)

if(any(as.vector(as.matrix(dat))<0)) doCCA <- FALSE else doCCA <- TRUE

doCCA <- FALSE

if(doFilt==TRUE){
	mycut <- min(abs(dat[which(dat!=0)]))*myfac
	good <- doFiltering(dat, cutoff=mycut, group=GROUPS, mypct=mypct)
	dat <- dat[good, ]
	mydata <- as.data.frame(dat)
	mydata$ID <- rownames(mydata)
	mydata <- mydata[, c("ID", colnames(dat))]
	write.table(mydata, file=paste("filtered.rel.freqs", counts.file, sep="."), sep="\t", row.names=FALSE, quote=FALSE)
	cond <- apply(dat, MARGIN=2, sum)
	dat <- dat[, which(cond>0)]
	meta <- meta[which(cond>0), ]; GROUPS <- GROUPS[which(cond>0)]
}

if(length(mylevels)>2){
	if(doCCA==TRUE) correspondence.analysis(MATRIX=dat, GROUPS=GROUPS, FILE=paste("adonis.CCA", "all.groups", study, "pdf", sep="."))
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
	
	if(doCCA==TRUE) correspondence.analysis(MATRIX=dat[, ind], GROUPS=GROUPS[ind], FILE=paste("adonis.CCA", tag, study, "pdf", sep="."))

	temdat <- t(as.matrix(dat[, ind]))
	rownames(temdat) <- names(DESCRIPTION)
	mytb <- differential.tests(MATRIX=temdat, GROUPS=GROUPS[ind], GROUP1=g1, GROUP2=g2, TEST="wilcox.test",
						PAIRED=mypaired, CUTOFF=2, COLNAME="feature", FILE=paste("wilcox.test", tag, study, sep="."), DESCRIPTION=DESCRIPTION)

	myfeatures <- rownames(dat)
	tem <- as.data.frame(temdat); tem$ID <- rownames(tem); tem <- tem[, c("ID", setdiff(colnames(tem), "ID"))]
	write.table(tem, file=paste("significant.rel.freqs", tag, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE) 
	
	if(mypaired==TRUE){
		pdf(file=paste("wilcox.test.paired.with.lines", tag, study, "pdf", sep="."))
		for(myfeature in myfeatures){
			plotPaired(dat=dat[, ind], feature=myfeature, group=GROUPS[ind],
					group.x=g1, group.y=g2, individual.factor=meta$donor[ind], DESCRIPTION=DESCRIPTION)
		}
		dev.off()
	}

	print(tag)
}

##
# si estamos con genes (TIGRFAM o KEGG), podemos poner la información de categ como una columna adicional en las tablas de resultados:
# (cargamos tabla que linca gene IDs con categorías o pathways)
categ <- NULL

if(length(grep("TIGR", rownames(dat)))==nrow(dat)){
	categ.tb <- read.table(file="tigrfam.roles.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
	categ <- categ.tb[, 3]
	names(categ) <- categ.tb[, 1]
}
if(length(grep("K", rownames(dat)))==nrow(dat)){
	categ.tb <- read.table(file="kegg.pathways.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
	categ <- categ.tb[, 3]
	names(categ) <- categ.tb[, 1]
}

if(!is.null(categ)){
	# leemos tablas a modificar (recien generadas)
	myfiles <- system("ls *.test.*.tsv", intern=TRUE)
	# asumimos que la columna con gene IDS es la primera
	for(file in myfiles){
		mytb <- read.table(file=file, header=TRUE, stringsAsFactors=FALSE, sep="\t", quote="")
		mytb$categ <- categ[mytb[, 1]]
		write.table(mytb, file=file, sep="\t", row.names=FALSE, quote=FALSE)
	}
}


