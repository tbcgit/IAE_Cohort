#!/usr/bin/env Rscript
###############

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

mylevels <- unlist(strsplit(mylevels, split=","))

ind <- which(is.element(GROUPS, mylevels))

if(do.arrow) myarrow <- TRUE else myarrow <- FALSE
principal.coordinates.analysis(MATRIX=as.dist(mydist[ind, ind]), GROUPS=GROUPS[ind], TITLE="", ELLIPSE=ellipse, FILE=paste(study, "pdf", sep="."),
	plot.arrows=myarrow, sum2lims=sum2lims)

