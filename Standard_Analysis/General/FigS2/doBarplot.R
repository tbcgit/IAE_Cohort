#!/usr/bin/env Rscript
###############

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
# myfactor ha de ser el nombre del factor conteniendo los grupos a agregar o el nombre de la columna con los sample names
myfactor <- as.character(parameters[3])
# si no aparece la casilla "other" en la legenda hay que aumentar este valor
mincutoff <- as.character(parameters[4])
tax.colum <- as.character(parameters[5])

if(mincutoff=="") mincutoff <- 0.6 else mincutoff <- as.numeric(mincutoff)

source("myfunctions.downstream.R")

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(dat) <- dat[, tax.colum]

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]
mysamples <- intersect(rownames(meta), colnames(dat))
mysamples <- mysamples[order(mysamples)]
mysamples <- c(mysamples[grep("I", mysamples)], mysamples[grep("A", mysamples)], mysamples[grep("E", mysamples)])
dat <- as.matrix(dat[, mysamples])
meta <- meta[mysamples, ]

if(all(meta[, myfactor]==rownames(meta))) byGroup <- FALSE else byGroup <- TRUE

GROUPS <- as.character(meta[, myfactor])
names(GROUPS) <- rownames(meta)

m <- dat
m <- 100*prop.table(m, margin=2)
m <- round(m*100)

# pasamos m a un data.frame cualitativo para poder aplicar las funciones de barplot de taxonomía
# con columnas "read", "feature", "sample"
func.tax <- c()
for(mysamp in colnames(m)){
	v <- c()
	for(mytax in rownames(m)){
		n <- m[mytax, mysamp]
		if(n>0) v <- c(v, rep(mytax, n))
	}
	temp <- data.frame(read=paste("r", 1:length(v)), feature=v, sample=rep(mysamp, length(v)))
	func.tax <- rbind(func.tax, temp)
}
func.tax$read <- as.character(func.tax$read)
func.tax$feature <- as.character(func.tax$feature)
func.tax$sample <- as.character(func.tax$sample)
func.tax$sample2 <- GROUPS[func.tax$sample]

func.tax$feature <- gsub("__", ".", func.tax$feature)
func.tax$feature <- gsub("_", ".", func.tax$feature)
func.tax$feature <- gsub("[(]", ".", func.tax$feature)
func.tax$feature <- gsub("[)]", ".", func.tax$feature)
func.tax$feature <- gsub("[[]", ".", func.tax$feature)
func.tax$feature <- gsub("[]]", ".", func.tax$feature)

if(byGroup==TRUE){
for(mygroup in unique(GROUPS)){
mysamples <- unique(func.tax$sample[which(func.tax$sample2==mygroup)])
names(mysamples) <- mysamples
myfilename <- paste("barplot", mygroup, sep="."); myfilename <- gsub(" ", ".", myfilename)
doSummary(mytb=func.tax, mincutoff=mincutoff, sampleColumn="sample", groups2plot=mysamples,
	GROUPS=mysamples, samplesORD=mysamples[length(mysamples):1], fixFileName=myfilename, mytaxlevel="genus")
print(mygroup)
}
} else {
	mysamples <- colnames(m)
	names(mysamples) <- mysamples
	myfilename <- "barplot.all.samples"
	doSummary(mytb=func.tax, mincutoff=mincutoff, sampleColumn="sample", groups2plot=mysamples,
	GROUPS=mysamples, samplesORD=mysamples[length(mysamples):1], fixFileName=myfilename, mytaxlevel="genus")
}

if(byGroup==TRUE){
	## agregando por grupos y calculando promedios de porcentajes
	mygroups <- unique(GROUPS)
	files <- gsub(" ", ".", paste("summary.taxa.barplot", mygroups, "csv", sep="."))
	mytaxa <- c()
	myphylum.v <- c()
	for(file in files){
		temp <- read.csv2(file=file, quote="", row.names=1)
		ind <- which(!is.element(rownames(temp), mytaxa))
		mytaxa <- c(mytaxa, rownames(temp)[ind])
		phfile <- sub("summary.taxa", "phylum", file)
		phfile <- sub("csv", "RData", phfile)
		load(file=phfile)
		myphylum.v <- c(myphylum.v, myphylum[ind])
		rm("myphylum")
	}
	ta.dfm.pct <- matrix(rep(0, length(mytaxa)*length(mygroups)), nrow=length(mygroups), ncol=length(mytaxa))
	rownames(ta.dfm.pct) <- mygroups; colnames(ta.dfm.pct) <- mytaxa
	for(group in mygroups){
		temp <- read.csv2(file=gsub(" ", ".", paste("summary.taxa.barplot", group, "csv", sep=".")), quote="", row.names=1)
		temp.means <- apply(temp, MARGIN=1, mean)
		ta.dfm.pct[group, names(temp.means)] <- temp.means
	}
	myfilename <- paste("summary.taxa", "by.group", myfactor, "csv", sep="."); myfilename <- gsub(" ", ".", myfilename)
	write.csv2(t(ta.dfm.pct), file=myfilename, quote=FALSE)
	pdf(file=sub("csv", "pdf", myfilename))
	if(all(is.element(c("I", "A", "E"), rownames(ta.dfm.pct)))){
		mybarplot(myphylum=myphylum.v, ta.dfm.pct=ta.dfm.pct[c("I", "A", "E"), ], "", mybottom=5)
	} else {
		mybarplot(myphylum=myphylum.v, ta.dfm.pct=ta.dfm.pct, "", mybottom=5)
	}
	dev.off()
}



