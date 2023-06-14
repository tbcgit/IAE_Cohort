#################
## Construcción de matrices de contribución taxonómica a KOs a partir del input de sqeezemeta 13.SQM_postselec_all_def.orftable.gz
# Antes seleccionamos columnas de interés para que ocupe menos el archivo que es muy muy grande (19GB).
# Nos quedamos solo con keggID, Tax y la parte de las muestras "raw_read_count":
# cut -f 9,10,487-721 13.SQM_postselec_all_def.orftable > mytable.tsv

dat <- read.table(file="mytable.tsv", sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
colnames(dat) <- sub("Raw read count ", "", colnames(dat))
kos <- scan(file="lista3.KOs.txt", what="character", sep="\n")

#######
# Lo hacemos por separado para las muestras del MTG y las del MTT
meta.file <- "metaMTT.tsv"
tag <- "RNA.species.contribution"

meta <- read.table(file=meta.file, sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)

mysamples <- intersect(meta$sample, colnames(dat))

mydat <- dat[, c("Tax", "KEGG ID", mysamples)]

ind2rm <- union(which(mydat[, "Tax"]==""), which(mydat[, "KEGG ID"]==""))
if(length(ind2rm)>0) mydat <- mydat[-ind2rm, ]

# Construimos todos estos niveles taxonómicos k_;p_;c_;o_;f_;g_;s_

tax.l <- strsplit(mydat$Tax, split=";")

tax <- unlist(lapply(tax.l, FUN=function(x){
			my.k <- x[1]
			ind <- grep("p_", x)
			if(length(ind)>0) my.p <- x[ind] else my.p <- "p_NA"
			ind <- grep("c_", x)
			if(length(ind)>0) my.c <- x[ind] else my.c <- "c_NA"
			ind <- grep("o_", x)
			if(length(ind)>0) my.o <- x[ind] else my.o <- "o_NA"
			ind <- grep("f_", x)
			if(length(ind)>0) my.f <- x[ind] else my.f <- "f_NA"
			ind <- grep("g_", x)
			if(length(ind)>0) my.g <- x[ind] else my.g <- "g_NA"
			ind <- grep("s_", x)
			if(length(ind)>0) my.s <- x[ind] else my.s <- "s_NA"
			return(paste(my.k, my.p, my.c, my.o, my.f, my.g, my.s, sep=";"))
		}))

mydat$Tax <- tax

for(ko in kos){
	#ind <- which(mydat[, "KEGG ID"]==ko)
	ind <- grep(ko, mydat[, "KEGG ID"])
	if(length(ind)>0){
		tem <- mydat[ind, ]
		cond <- which(apply(tem[, 3:ncol(tem)], MARGIN=1, sum)>0)
		if(length(cond)>0){
			tem <- tem[cond, ]
			myrownames <- unique(tem$Tax)
			myncol <- length(mysamples)
			mynrow <- length(myrownames)
			M <- matrix(rep(0, myncol*mynrow), ncol=myncol, nrow=mynrow)
			rownames(M) <- myrownames; colnames(M) <- mysamples
			for(tax in myrownames) M[tax, ] <- apply(tem[which(tem$Tax==tax), colnames(M)], MARGIN=2, sum)
			if(!is.null(nrow(M))) M <- as.data.frame(M) else M <- as.data.frame(t(as.matrix(M)))
			M$ID <- rownames(M)
			M <- M[, c("ID", mysamples)]
			write.table(M, file=paste(ko, tag, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)
		}
	}
	print(ko)
}

