#!/usr/bin/env Rscript
###############
#./getPositiveMatrix.R ANCOMBC.normalized.counts.genus.tsv metadata.unique.tsv ANCOMBC.norm.pos.counts.genus.tsv ID

parameters <- commandArgs(trailingOnly=TRUE)

input.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
out.file <- as.character(parameters[3])
descr.column <- as.character(parameters[4])

dat <- read.table(file=input.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
DESCRIPTION <- dat[, descr.column]
rownames(dat) <- dat[, 1]

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, row.names=1, check.names=FALSE)

dat <- dat[, intersect(rownames(meta), colnames(dat))]

if(any(dat<0)) dat <- dat + abs(min(dat, na.rm=TRUE))

dat <- as.data.frame(dat)

dat$ID <- rownames(dat)
dat$DESCRIPTION <- DESCRIPTION
dat <- dat[, c("ID", setdiff(colnames(dat), c("ID", "DESCRIPTION")), "DESCRIPTION")]

write.table(dat, file=out.file, sep="\t", row.names=FALSE, quote=FALSE)


