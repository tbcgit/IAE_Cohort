#!/usr/bin/env Rscript
###############

parameters <- commandArgs(trailingOnly=TRUE)

counts.file <- as.character(parameters[1])
meta.file <- as.character(parameters[2])
myfactor <- as.character(parameters[3])
study <- as.character(parameters[4])
columns.char <- as.character(parameters[5])

columns.char <- unlist(strsplit(columns.char, split=","))
descr.column <- columns.char[1]

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
DESCRIPTION <- dat[, descr.column]; names(DESCRIPTION) <- rownames(dat) <- dat[, 1]

if(length(columns.char)>1) dat.extra <- dat[, columns.char[2:length(columns.char)]]

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
rownames(meta) <- meta[, 1]
mysamples <- intersect(rownames(meta), colnames(dat))
dat <- as.matrix(dat[, mysamples])
meta <- meta[mysamples, ]

GROUPS <- as.character(meta[, myfactor])
names(GROUPS) <- mysamples

mygroups <- unique(GROUPS)
M <- matrix(rep(0, length(mygroups)*nrow(dat)), nrow=nrow(dat), ncol=length(mygroups))
rownames(M) <- rownames(dat)
colnames(M) <- mygroups

# normalizamos
dat <- prop.table(dat, margin=2)*100
# pasamos a counts otra vez
for(i in 1:ncol(dat)) dat[, i] <- round(dat[, i]/min(dat[which(dat[, i]!=0), i]))

for(mygroup in mygroups){
	ind <- which(GROUPS==mygroup)
	if(length(ind)>1) M[, mygroup] <- round(apply(dat[, ind], MARGIN=1, mean)) else M[, mygroup] <- dat[, ind]
}

M <- as.data.frame(M)

M$ID <- rownames(M); M$DESCRIPTION <- DESCRIPTION[rownames(M)]

M <- M[, c("ID", mygroups, "DESCRIPTION")]

if(length(columns.char)>1) M <- cbind(M, dat.extra)

meta <- meta[which(!duplicated(GROUPS)), ]
meta[, 1] <- meta[, myfactor]

write.table(M, file=paste("mean.counts", study, "tsv", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
write.table(meta, file=paste("metadata", study, "tsv", sep="."), sep="\t", quote=FALSE, row.names=FALSE)

