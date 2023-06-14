#!/usr/bin/env Rscript
###############

parameters <- commandArgs(trailingOnly=TRUE)

otu.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
out.file <- as.character(parameters[3])

library(vegan)

dat <- read.table(file=otu.file, sep="\t", header=TRUE, quote="", row.names=1, stringsAsFactors=FALSE, check.names=FALSE)

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, row.names=1, check.names=FALSE)

dat <- t(dat[, intersect(rownames(meta), colnames(dat))])

bray.dist <- vegdist(x=dat, method="euclidean")

write.csv(as.matrix(bray.dist), file=out.file)


