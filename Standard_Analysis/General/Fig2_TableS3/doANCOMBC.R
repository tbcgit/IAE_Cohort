#!/usr/bin/env Rscript
#/scratch/software/Rstats/bin/Rscript
###############
# Este script solo sirve para comparar 2 grupos.
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
# El penúltimo parámetro, si no es un estudio pareado, ha de valer "group". Si es un estudio pareado ha de valer "group + donor".
#
# ./doANCOMBC.R counts.asv.tsv meta.asv.tsv "GVHDs,week" "no.pre,yes.pre" "sinfilt" FALSE 0 0 TAXONOMY "group" FALSE

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

if(length(grep("donor", myformula))>0) mypaired <- TRUE else mypaired <- FALSE

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

if(mypaired==TRUE) meta$donor <- as.character(meta$donor)

dat.pct <- prop.table(dat, margin=2)*100

dat.nofilt <- dat

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
} else { good <- rownames(dat) }

g1 <- mylevels[1]
g2 <- mylevels[2]
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

meta.ancom <- meta[ind, ]
meta.ancom$group <- GROUPS[ind]
OTU <- otu_table(dat.nofilt[, ind], taxa_are_rows=TRUE)
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

# filtramos después de transformar
dat.trans <- dat.trans[good, ]

if(structzero==TRUE){
	structural_zeros <- as.data.frame(out$zero_ind)
	structural_zeros$DESCRIPTION <- DESCRIPTION[rownames(structural_zeros)]
	structural_zeros$ID <- rownames(structural_zeros)
	structural_zeros <- structural_zeros[, c("ID", setdiff(colnames(structural_zeros), c("ID", "DESCRIPTION")), "DESCRIPTION")]
	cond1 <- apply(structural_zeros[, c(2,3)], MARGIN=1, FUN=function(x) if(all(x==TRUE)) return(TRUE) else return(FALSE))
	cond2 <- apply(structural_zeros[, c(2,3)], MARGIN=1, FUN=function(x) if(any(x==TRUE)) return(TRUE) else return(FALSE))
	structural_zeros <- rbind(structural_zeros[which(cond1==TRUE), ], structural_zeros[which(cond2==TRUE),], structural_zeros[which(cond2==FALSE),])
	write.table(structural_zeros, file=paste("structural_zeros", tag, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
}

df <- data.frame(feature=rownames(out$feature_table), ANCOMBC.pval=res$p_val[, 1], ANCOMBC.adjpval.fdr=res$q_val[, 1], stringsAsFactors=FALSE)
df$DESCRIPTION <- DESCRIPTION[df$feature]
if(structzero==TRUE) df <- cbind(df, structural_zeros[df$feature, c(2,3)])

mytb <- differential.tests(MATRIX=dat.pct[, ind], GROUPS=GROUPS[ind], GROUP1=g1, GROUP2=g2, TEST="wilcox.test",
						PAIRED=mypaired, CUTOFF=0.1, COLNAME="feature", FILE=paste("ANCOMBC.test", tag, study, sep="."), DESCRIPTION=DESCRIPTION)
mytb$feature <- as.character(mytb$feature)
rownames(mytb) <- mytb$feature
mytb <- mytb[, c("mean1","mean2","log2FC")]

# Filtramos tabla de p-valores de ANCOMBC y recalculamos p-valor ajustado
rownames(df) <- df$feature
df <- df[good, ]
df[, "ANCOMBC.adjpval.fdr"] <- p.adjust(df[, "ANCOMBC.pval"], method="fdr")
myord <- order(df$ANCOMBC.pval, decreasing=FALSE)
df <- df[myord, ]

mytb <- mytb[rownames(df), ]
df <- cbind(df, mytb)
df <- df[, c(setdiff(colnames(df), "DESCRIPTION"), "DESCRIPTION")]

myfeatures <- as.character(df[which(df$ANCOMBC.pval<=0.1), "feature"])

# pintamos boxplots con datos transformados (machacaremos la gráfica de differential.tests)
G1 <- names(GROUPS)[which(GROUPS==g1)]
G2 <- names(GROUPS)[which(GROUPS==g2)]
items.group1 <- which(is.element(colnames(dat.trans), G1))
items.group2 <- which(is.element(colnames(dat.trans), G2))

pdf(file=paste("ANCOMBC.test", tag, study, "pdf", sep="."))
for(sign.row in myfeatures){
    	values <- list(dat.trans[sign.row, items.group1], dat.trans[sign.row, items.group2])
    	names(values) <- c(g1, g2)
    	mymain <- DESCRIPTION[sign.row]
    	beeswarm(values, pch=16, col=c("red","blue"), las=1, xlab="Condition", ylab=sign.row, main=mymain, cex.lab=0.6, cex.main=0.5)
    	graphics::boxplot(values, las=1,  add=TRUE, col=NULL)
    	legend("topright", legend=c(paste("ANCOMBC.pval:"         , format(df[sign.row, "ANCOMBC.pval"    ], digits=2)),
                                paste("ANCOMBC.adjpval.fdr:", format(df[sign.row, "ANCOMBC.adjpval.fdr"], digits=2))), cex=0.7)
}
dev.off()
	
if(mypaired==TRUE){
	pdf(file=paste("ANCOMBC.test.paired.with.lines", tag, study, "pdf", sep="."))
	for(myfeature in myfeatures){
		plotPaired(dat=dat.trans, feature=myfeature, group=GROUPS[colnames(dat.trans)],
				group.x=g1, group.y=g2, individual.factor=meta[colnames(dat.trans), "donor"],
				DESCRIPTION=DESCRIPTION, p.val=df[myfeature, "ANCOMBC.pval"])
	}
	dev.off()
}

df$feature <- origIDs[df$feature]
write.table(df, file=paste("ANCOMBC.test", tag, study, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)

dat.trans <- as.data.frame(dat.trans)
dat.trans$ID <- origIDs[rownames(dat.trans)]
dat.trans <- dat.trans[, c("ID", setdiff(colnames(dat.trans), "ID"))]
dat.trans$DESCRIPTION <- DESCRIPTION[rownames(dat.trans)]
write.table(dat.trans, file=paste("ANCOMBC.normalized", counts.file, sep="."), sep="\t", row.names=FALSE, quote=FALSE)

# generamos input boruta
tem <- data.frame(feature=df$feature, p.value=rep(0, nrow(df)), adj.p.value=rep(0, nrow(df)), stringsAsFactors=FALSE)
write.table(tem, file=paste("input.boruta", counts.file, sep="-"), sep="\t", row.names=FALSE, quote=FALSE)

print(tag)

