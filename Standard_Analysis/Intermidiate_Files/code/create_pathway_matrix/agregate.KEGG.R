################
## Agrega filas de matriz abundancias de KEGG models a los niveles C y B usando la tabla
## ko00001.keg.tsv (previamente generada con script.KEGG.R)

## devuelve matrices de abundancia a los niveles D, C y B, con primera columna ID, e incluyendo columna DESCRIPTION

kegg.matrix.file <- "SQM_postselec_all_def.KO.tpm.tsv"
categ.file <- "ko00001.keg.tsv"
out.name <- "tpm"

kegg.matrix <- read.table(file=kegg.matrix.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE, row.names=1)

categ <- read.table(file=categ.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE, comment.char="")
categ$C.id <- as.character(categ$C.id); categ$B.id <- as.character(categ$B.id)

tem <- categ[which(!duplicated(categ$D.id)), ]
D <- tem$D
names(D) <- tem$D.id

tem <- categ[which(!duplicated(categ$C.id)), ]
C <- tem$C
names(C) <- tem$C.id

tem <- categ[which(!duplicated(categ$B.id)), ]
B <- tem$B
names(B) <- tem$B.id

# D
D.mat <- kegg.matrix
D.mat[, "ID"] <- rownames(D.mat)
D.mat[, "DESCRIPTION"] <- D[D.mat[, "ID"]]
mysamples <- setdiff(colnames(D.mat), c("ID", "DESCRIPTION"))
D.mat <- D.mat[, c("ID", mysamples, "DESCRIPTION")]

# C
C.mat <- matrix(rep(0, length(C)*length(mysamples)), nrow=length(C), ncol=length(mysamples))
rownames(C.mat) <- names(C); colnames(C.mat) <- mysamples
for(i in 1:nrow(C.mat)){
	mylab <- rownames(C.mat)[i]
	mykeggs <- intersect(unique(categ$D.id[which(categ$C.id==mylab)]), rownames(D.mat))
	if(length(mykeggs)>0) C.mat[mylab, mysamples] <- apply(D.mat[mykeggs, mysamples], MARGIN=2, sum)
}
cond <- apply(C.mat, MARGIN=1, sum)
C.mat <- C.mat[which(cond>0), ]
C.mat <- as.data.frame(C.mat)
C.mat[, "ID"] <- rownames(C.mat)
C.mat[, "DESCRIPTION"] <- C[C.mat[, "ID"]]
C.mat <- C.mat[, c("ID", mysamples, "DESCRIPTION")]

# B
B.mat <- matrix(rep(0, length(B)*length(mysamples)), nrow=length(B), ncol=length(mysamples))
rownames(B.mat) <- names(B); colnames(B.mat) <- mysamples
for(i in 1:nrow(B.mat)){
	mylab <- rownames(B.mat)[i]
	mykeggs <- intersect(unique(categ$D.id[which(categ$B.id==mylab)]), rownames(D.mat))
	if(length(mykeggs)>0) B.mat[mylab, mysamples] <- apply(D.mat[mykeggs, mysamples], MARGIN=2, sum)
}
cond <- apply(B.mat, MARGIN=1, sum)
B.mat <- B.mat[which(cond>0), ]
B.mat <- as.data.frame(B.mat)
B.mat[, "ID"] <- rownames(B.mat)
B.mat[, "DESCRIPTION"] <- B[B.mat[, "ID"]]
B.mat <- B.mat[, c("ID", mysamples, "DESCRIPTION")]

write.table(D.mat, file=paste("kegg", out.name, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
write.table(C.mat, file=paste("pathway.kegg", out.name, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
write.table(B.mat, file=paste("subcategory.kegg", out.name, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)


