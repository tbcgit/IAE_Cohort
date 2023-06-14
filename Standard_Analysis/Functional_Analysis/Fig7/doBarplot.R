#!/usr/bin/env Rscript
###############
##### doBarplot
# Se aplica sobre la matriz a nivel de "family", "genus" o "species" que se obtiene con splitTaxLevels.
# (parámetro 3) ha de ser el nombre del factor conteniendo los grupos a agregar, o el nombre de la columna en metadata con los sample names, si lo
#			que queremos es representar todas las muestras del metadata (y en el orden dado en el metadata) en el mismo plot.
# Si (parámetro 3) es el nombre de un factor con los grupos a agregar se genera:
#	(1) un plot por grupo con todas las muestras del grupo representadas en cada plot
#	(2) un plot representando cada grupo en una única barra con sus muestras agregadas
# NOTA: El agregamiento para cada taxón representado consiste en calcular el promedio de los porcentajes de ese taxón en todas las muestras del grupo.
# (parámetro 4) Si no aparece la casilla "other" en la legenda se puede aumentar el valor de este parámetro.
# (parámetro 6) Ha de ser uno de entre "family", "genus" o "species"

#./doBarplot.R g_counts.tsv Metadata_without_neg.tsv "Group" 0.6 "Taxon" "genus"

#./doBarplot.R g_counts.tsv Metadata_without_neg.tsv "Samples" 0.6 "Taxon" "genus"


parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
# myfactor ha de ser el nombre del factor conteniendo los grupos a agregar o el nombre de la columna con los sample names
myfactor <- as.character(parameters[3])
# si no aparece la casilla "other" en la legenda hay que aumentar este valor
mincutoff <- as.character(parameters[4])
tax.colum <- as.character(parameters[5])
mytaxlevel <- as.character(parameters[6])

if(mincutoff=="") mincutoff <- 0.6 else mincutoff <- as.numeric(mincutoff)

source("myfunctions.downstream.R")

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(dat) <- dat[, tax.colum]

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]
mysamples <- intersect(rownames(meta), colnames(dat))
dat <- as.matrix(dat[, mysamples])
meta <- meta[mysamples, ]

if(all(meta[, myfactor]==rownames(meta))) byGroup <- FALSE else byGroup <- TRUE

GROUPS <- as.character(meta[, myfactor])
names(GROUPS) <- rownames(meta)

m <- dat
m <- 100*prop.table(m, margin=2)
m <- round(m*100)

########### *** arreglo ***
## Los números que se pegan al final de la taxonomía repetida se corresponden con el número de fila en la tabla de entrada.
tem.l <- strsplit(rownames(m), split=";")
len <- length(tem.l[[1]])
if(mytaxlevel=="species") tem <- unlist(lapply(tem.l, FUN=function(x) paste(x[4:6], collapse="|")))
if(mytaxlevel=="genus") tem <- unlist(lapply(tem.l, FUN=function(x) paste(x[4:5], collapse="|")))
if(mytaxlevel=="family") tem <- unlist(lapply(tem.l, FUN=function(x) x[4]))
inds <- which(duplicated(tem))
if(length(inds)>0){
	myinds <- c()
	for(ind in inds) myinds <- c(myinds, which(tem==tem[ind]))
	myinds <- unique(myinds)
	for(ind in myinds){
		tem.l[[ind]][4] <- paste(tem.l[[ind]][2], tem.l[[ind]][4], sep="|")
		tem.l[[ind]][len] <- paste(tem.l[[ind]][len], ind, sep="_")
	}
	rownames(m) <- unlist(lapply(tem.l, FUN=function(x) paste(x, collapse=";")))
}
###########

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
	GROUPS=mysamples, samplesORD=mysamples[length(mysamples):1], fixFileName=myfilename, mytaxlevel=mytaxlevel)
print(mygroup)
}
} else {
	mysamples <- colnames(m)
	names(mysamples) <- mysamples
	myfilename <- "barplot.all.samples"
	doSummary(mytb=func.tax, mincutoff=mincutoff, sampleColumn="sample", groups2plot=mysamples,
	GROUPS=mysamples, samplesORD=mysamples[length(mysamples):1], fixFileName=myfilename, mytaxlevel=mytaxlevel)
}

if(byGroup==TRUE){
	## agregando por grupos y calculando promedios de porcentajes
	mygroups <- unique(GROUPS)
	files <- gsub(" ", ".", paste("summary.taxa.barplot", mygroups, "csv", sep="."))
	mytaxa <- c()
	myphylum.v <- c()
	for(file in files){
		temp <- read.csv(file=file, quote="", row.names=1)
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
		temp <- read.csv(file=gsub(" ", ".", paste("summary.taxa.barplot", group, "csv", sep=".")), quote="", row.names=1)
		temp.means <- apply(temp, MARGIN=1, mean)
		ta.dfm.pct[group, names(temp.means)] <- temp.means
	}
	myfilename <- paste("summary.taxa", "by.group", myfactor, "csv", sep="."); myfilename <- gsub(" ", ".", myfilename)
	write.csv(t(ta.dfm.pct), file=myfilename, quote=FALSE)
	pdf(file=sub("csv", "pdf", myfilename))
	mybarplot(myphylum=myphylum.v, ta.dfm.pct=ta.dfm.pct, "", mybottom=5)
	dev.off()
}



