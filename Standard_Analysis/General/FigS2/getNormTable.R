#!/usr/bin/env Rscript
###############
# ./getNormTable.R averaged.counts.MTGs.tsv averaged.metadata.tsv averaged.norm.MTGs.tsv 1 "ID,DESCRIPTION" FALSE

parameters <- commandArgs(trailingOnly=TRUE)

otu.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
out.file <- as.character(parameters[3])
type <- as.numeric(as.character(parameters[4]))
columns.char <- as.character(parameters[5])
doLog <- as.logical(as.character(parameters[6]))

columns.char <- unlist(strsplit(columns.char, split=","))

dat <- read.table(file=otu.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)

if(nrow(dat)>1){
rownames(dat) <- dat[, 1]

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, row.names=1, check.names=FALSE)

mysamples <- intersect(rownames(meta), colnames(dat))

dat.M <- dat[, mysamples]

if(type==1){
	# TSS
	mat <- prop.table(as.matrix(dat.M), margin=2)*100
	if(doLog){
		mat <- mat+1
		mat <- log(mat)
	}
}

if(type==2){
	# CSS
	library(metagenomeSeq)
	counts <- dat
	counts$OTU <- dat[, 1]
	counts <- counts[, c("OTU", mysamples)]
	write.table(counts, file="counts.tsv", quote=FALSE, sep="\t", row.names=FALSE)
	# Loading count data
	mydat <- load_meta("counts.tsv")
	# Creating a MRexperiment object
	obj <- newMRexperiment(mydat$counts)
	## Normalisation
	# determine the proper percentile to normalize
	p <- cumNormStat(obj)
	# calculate the scaling factors
	obj <- cumNorm(obj, p=p)
	# call normalized counts
	mat <- MRcounts(obj, norm=TRUE, log=TRUE)
}

if(type==3){
	# centered log-ratio (CLR) transformation
	mat <- as.matrix(dat.M)
	mat <- prop.table(mat, margin=2)*100
	# Cambiamos la escala para que el máximo valor sea 100, de otro modo, si los valores son muy grandes y hay muchas muestras
	# no se puede calcular el producto de toda una fila, R se queda sin precisión y da infinito.
	#mat <- 1 + 100*mat/max(mat)
	#for(i in 1:ncol(mat)){
	#	G <- prod(mat[, i])**(1/nrow(mat))
	#	mat[, i] <- log(mat[, i]/G)
	#}
	library(compositions)
	mat <- as.data.frame(t(clr(t(mat))))
	
}

if(type!=1&type!=2&type!=3){
	print("Tienes que dar al parámetro 4 el valor 1, 2 ó 3")
	mat <- dat.M
}

mat <- as.data.frame(mat)

if(length(columns.char)>1){
	mat[, columns.char] <- dat[, columns.char]
	mat <- mat[, c(columns.char[1], mysamples, columns.char[2:length(columns.char)])]
} else {
	mat[, columns.char] <- dat[, columns.char]
	mat <- mat[, c(columns.char[1], mysamples)]
}

write.table(mat, file=out.file, sep="\t", row.names=FALSE, quote=FALSE)
}

