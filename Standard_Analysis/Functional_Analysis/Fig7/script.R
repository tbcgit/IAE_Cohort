#############
##

datDNA <- read.table(file="agregado.g_K01696.DNAcontribution.tsv", sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
datDNA <- as.matrix(datDNA[, c("A", "E", "I")])

datRNA <- read.table(file="agregado.g_K01696.RNAcontribution.tsv", sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
datRNA <- as.matrix(datRNA[, c("A", "E", "I")])

myrownames <- union(rownames(datDNA), rownames(datRNA))
mycolnames <- c(paste(c("A", "E", "I"), "DNA", sep="."), paste(c("A", "E", "I"), "RNA", sep="."))
M <- matrix(rep(0, length(myrownames)*length(mycolnames)), nrow=length(myrownames), ncol=length(mycolnames))
rownames(M) <- myrownames
colnames(M) <- mycolnames

M[rownames(datDNA), c("A.DNA", "E.DNA", "I.DNA")] <- datDNA[rownames(datDNA), c("A", "E", "I")]
M[rownames(datRNA), c("A.RNA", "E.RNA", "I.RNA")] <- datRNA[rownames(datRNA), c("A", "E", "I")]

M <- as.data.frame(M)
M$ID <- rownames(M)
M <- M[, c("ID", setdiff(colnames(M), "ID"))]
write.table(M, file="g_K01696.tsv", sep="\t", row.names=FALSE, quote=FALSE)

#####

# TnaA
dat1667 <- read.table(file="g_K01667.tsv", sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
dat1667 <- as.matrix(dat1667)
colnames(dat1667) <- paste(colnames(dat1667), "TnaA", sep=".")

# TrpB
dat1696 <- read.table(file="g_K01696.tsv", sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
dat1696 <- as.matrix(dat1696)
colnames(dat1696) <- paste(colnames(dat1696), "TrpB", sep=".")

myrownames <- union(rownames(dat1667), rownames(dat1696))
mycolnames <- c(colnames(dat1667), colnames(dat1696))
M <- matrix(rep(0, length(myrownames)*length(mycolnames)), nrow=length(myrownames), ncol=length(mycolnames))
rownames(M) <- myrownames
colnames(M) <- mycolnames

M[rownames(dat1667), colnames(dat1667)] <- dat1667[rownames(dat1667), colnames(dat1667)]
M[rownames(dat1696), colnames(dat1696)] <- dat1696[rownames(dat1696), colnames(dat1696)]

M <- as.data.frame(M)
M$ID <- rownames(M)
M <- M[, c("ID", setdiff(colnames(M), "ID"))]
write.table(M, file="g_counts.tsv", sep="\t", row.names=FALSE, quote=FALSE)

