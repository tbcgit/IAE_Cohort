#!/usr/bin/env Rscript
###############
# Si mylevels = "" se entiende que myfactors (variable explicativa) se corresponde con una variable continua. En ese caso
# se pinta el pcoa fabricando 4 grupos a partir de la variable continua, y usando flechas desde el origen hasta las muestras. 
# ./doAdonisFromDistances.R genus.euclidean.dissimilarities.csv metadata_FC.txt "log2FC" "" FC FALSE TRUE 0

parameters <- commandArgs(trailingOnly=TRUE)

dist.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
myfactors <- as.character(parameters[3])
mylevels <- as.character(parameters[4])
study <- as.character(parameters[5])
ellipse <- as.logical(as.character(parameters[6]))
do.arrow <- as.logical(as.character(parameters[7]))
sum2lims <- as.numeric(parameters[8])

source("myfunctions.downstream.R")

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]

mydist <- read.csv(file=dist.file, row.names=1, check.names=FALSE)
mysamples <- intersect(rownames(meta), rownames(mydist))
mydist <- mydist[mysamples, mysamples]
meta <- meta[mysamples, ]

myfactors.l <- strsplit(myfactors, split=",")
myfactor1 <- myfactors.l[[1]][1]
if(length(myfactors.l[[1]])==2) myfactor2 <- myfactors.l[[1]][2] else myfactor2 <- ""
if(myfactor2!="") GROUPS <- paste(meta[, myfactor1], meta[, myfactor2], sep=".") else GROUPS <- as.character(meta[, myfactor1])
names(GROUPS) <- mysamples

if(mylevels!=""){
	mylevels <- unlist(strsplit(mylevels, split=","))
	ind <- which(is.element(GROUPS, mylevels))
} else {
	GROUPS <- as.numeric(GROUPS); names(GROUPS) <- mysamples
	ind <- 1:length(GROUPS)
}

if(do.arrow) myarrow <- TRUE else myarrow <- FALSE

myD <- as.dist(mydist[ind, ind])

GROUPS <- GROUPS[ind]

mydf <- data.frame(sample=names(GROUPS), GROUPS=GROUPS)

adonis.res <- adonis(myD ~ GROUPS, data=mydf, perm=600)

save(adonis.res, file=paste(study, "RData", sep="."))

adonis.p.value <- format(adonis.res$aov.tab["GROUPS", "Pr(>F)"], digits=3)

adonis.R2 <- format(adonis.res$aov.tab["GROUPS", "R2"], digits=3)

if(mylevels==""){
	GROUPS <- as.character(cut(x=GROUPS, breaks=4)); names(GROUPS) <- mysamples
}

principal.coordinates.analysis(MATRIX=myD, GROUPS=GROUPS, TITLE=paste("permanova p-value: ", adonis.p.value, "; R-squared: ", adonis.R2, sep=""),
	ELLIPSE=ellipse, FILE=paste(study, "pdf", sep="."), plot.arrows=myarrow, sum2lims=sum2lims)

