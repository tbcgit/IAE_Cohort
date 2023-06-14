#!/usr/bin/env Rscript
#/scratch/software/Rstats/bin/Rscript
###############
# Este script sirve para normalizar la matriz de abundancias según el método ANCOMBC cuando queremos que la matriz contenga más de 2 grupos.
#
# La matriz de entrada ha de ser una matriz de frecuencias absolutas.
#
# Si el último parámetro es TRUE, se identifican ceros estructurales (casi todo ceros en un grupo)
# y se les asigna por defecto un p-valor y p-valor ajustado de cero. Si es FALSE, se calcula el p-valor de esos structural zeros
# como si fueran cualquier otra variable, obteniéndose un valor muy cercano a cero pero no cero.
# Ellos en el paper intentan vender eso de que los structural zeros deberían ser automáticamente significativos (de ahí lo de p-valor cero),
# pero en la función de R el tratamiento por defecto es calcular su p-valor como si se tratara de cualquier otra variable.
# Parece como si eso de dar un p-valor de cero a los structural zeros fuera una decisión que asumieron en versiones anteriores
# de la herramienta y que ahora se hayan arrepentido, aunque sigan dando la opción de hacerlo por coherencia con sus métodos antiguos.
#
# El penúltimo parámetro ha de valer "group". Podríamos incluir covariables adicionales con la forma "group + cov1 + cov2 + ...".
#
# ./doANCOMBC_MT2G.R counts.species.tsv metadata.10K.updated.tsv "treatment" "BAL,BAS,BIO,CEP,ES,EX" "sinfilt" FALSE 0 0 TAXONOMY "group" FALSE

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
myformula <- as.character(parameters[10])
structzero <- as.logical(as.character(parameters[11]))

source("myfunctions.downstream.R")

library(ANCOMBC)
library(phyloseq)

dat <- read.table(file=counts.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
origIDs <- dat[, 1]
dat[, 1] <- gsub(";", ".", dat[, 1])
dat[, 1] <- gsub("_", ".", dat[, 1])
dat[, 1] <- gsub("-", ".", dat[, 1])
dat[, 1] <- gsub("[(]", ".", dat[, 1])
dat[, 1] <- gsub("[)]", ".", dat[, 1])
dat[, 1] <- gsub("[[]", ".", dat[, 1])
dat[, 1] <- gsub("[]]", ".", dat[, 1])
dat[, 1] <- gsub(" ", ".", dat[, 1])
dat[, 1] <- gsub(":", ".", dat[, 1])
dat[, 1] <- gsub("/", ".", dat[, 1])
dat[, 1] <- gsub("[|]", ".", dat[, 1])
names(origIDs) <- dat[, 1]

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

dat.pct <- prop.table(dat, margin=2)*100

if(doFilt==TRUE){
	mycut <- min(abs(dat.pct[which(dat!=0)]))*myfac
	good <- doFiltering(dat.pct, cutoff=mycut, group=GROUPS, mypct=mypct)
	dat <- dat[good, ]; dat.pct <- dat.pct[good, ]
	mydata <- as.data.frame(dat)
	mydata$ID <- rownames(mydata)
	mydata <- mydata[, c("ID", colnames(dat))]
	write.table(mydata, file=paste("filtered.rel.freqs", counts.file, sep="."), sep="\t", row.names=FALSE, quote=FALSE)
	cond <- apply(dat, MARGIN=2, sum)
	dat <- dat[, which(cond>0)]; dat.pct <- dat.pct[, which(cond>0)]
	meta <- meta[which(cond>0), ]; GROUPS <- GROUPS[which(cond>0)]
}

meta.ancom <- meta
meta.ancom$group <- GROUPS
OTU <- otu_table(dat, taxa_are_rows=TRUE)
META <- sample_data(meta.ancom)
physeq <- phyloseq(OTU, META)

#### Running ancombc function.
# conserve: whether to use a conservative variance estimate of the test statistic. It is recommended
# if the sample size is small and/or the number of differentially abundant taxa is believed to be large. Default is FALSE.
out <- ancombc(phyloseq=physeq, formula=myformula, 
              p_adj_method="fdr", zero_cut=0.95, lib_cut=0, 
              group="group", struc_zero=structzero, neg_lb=FALSE, tol=1e-5, 
              max_iter=100, conserve=FALSE, global=FALSE)

res <- out$res

#### Bias-adjusted abundances
samp_frac <- out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] <- 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn <- log(out$feature_table + 1) 
# Adjust the log observed abundances
dat.trans <- t(t(log_obs_abn) - samp_frac)

if(structzero==TRUE){
	structural_zeros <- as.data.frame(out$zero_ind)
	structural_zeros$DESCRIPTION <- DESCRIPTION[rownames(structural_zeros)]
	structural_zeros$ID <- rownames(structural_zeros)
	structural_zeros <- structural_zeros[, c("ID", setdiff(colnames(structural_zeros), c("ID", "DESCRIPTION")), "DESCRIPTION")]
	write.table(structural_zeros, file=paste("structural_zeros", study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
}

dat.trans <- as.data.frame(dat.trans)
dat.trans$ID <- origIDs[rownames(dat.trans)]
dat.trans <- dat.trans[, c("ID", setdiff(colnames(dat.trans), "ID"))]
dat.trans$DESCRIPTION <- DESCRIPTION[rownames(dat.trans)]
write.table(dat.trans, file=paste("ANCOMBC.normalized", counts.file, sep="."), sep="\t", row.names=FALSE, quote=FALSE)



