#####################
## Downstream para los outputs de complexCruncher

library("XLConnect")

taxColumn <- "GENUS"

########
## aplicar CCA y wilcox (Control vs LS) a la matriz [lineageX(RSI de todas las muestras)]

mygroup <- "group"

# el metadata ha de tener las mismas muestras que aparecen en el data set de complexCruncher (no necesario el mismo orden)
meta <- read.table(file="meta.tsv", sep="\t", header=TRUE, quote="", stringsAsFactors=FALSE, check.names=FALSE)

ind2rm <- which(is.na(meta[, mygroup]))
if(length(ind2rm)>0) meta <- meta[-ind2rm, ]

mysamples <- meta[, 1]
GROUPS <- meta[, mygroup]
names(GROUPS) <- mysamples

groups.l <- vector("list", length(unique(GROUPS)))
names(groups.l) <- unique(GROUPS)
for(g in unique(GROUPS)) groups.l[[g]] <- mysamples[which(GROUPS==g)]

# Si trabajamos sólo con la intersección de los linajes de todas las tablas al final obtenemos muy pocas,
# de modo que trabajaremos con la unión e imputaremos un linaje NA en una muestra con el valor de RSI neutral de 0.5.
# Esta tabla está en la carpeta corrank de resultados de complexCruncher.
wb <- loadWorkbook("counts.genus.individuals.hA_rank.xlsx")
myvars <- readWorksheet(wb, sheet = 1)[, taxColumn]
for(s in mysamples){
	df <- readWorksheet(wb, sheet = s)
	myvars <- c(myvars, df[, taxColumn])
	print(summary(df$RSI))
}
myvars <- unique(myvars)

df.l <- vector("list", length(mysamples))
names(df.l) <- mysamples
for(s in mysamples) df.l[[s]] <- readWorksheet(wb, sheet = s)

M <- matrix(rep(NA, length(myvars)*length(mysamples)), nrow=length(myvars), ncol=length(mysamples))
rownames(M) <- myvars
colnames(M) <- mysamples
for(myvar in myvars){
	for(mysample in mysamples){
		rownames(df.l[[mysample]]) <- df.l[[mysample]][, taxColumn]
		M[myvar, mysample] <- df.l[[mysample]][myvar, "RSI"]
	}
}

# imputamos
for(j in 1:ncol(M)){
	ind <- which(is.na(M[, j]))
	#m <- min(M[, j], na.rm=TRUE)
	#M[ind, j] <- runif(n=length(ind), min=m/2, max=m)
	m <- mean(M[, j], na.rm=TRUE)
	M[ind, j] <- runif(n=length(ind), min=m-m/5, max=m+m/5)
}

M <- as.data.frame(M)
M$ID <- rownames(M)
M <- M[, c("ID", setdiff(colnames(M), "ID"))]
write.table(M, file="RSI.matrix.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# Hacemos análisis desde la terminal linux:
# ./getEuclidianDist.R RSI.matrix.tsv meta.all.LS.tsv mydist.csv
# ./doPCOA.R mydist.csv meta.all.LS.tsv "Gender" "F,M" pcoa.Gender TRUE FALSE 0.15
# ./doTest.R RSI.matrix.tsv meta.all.LS.tsv "Gender" "F,M" Gender FALSE 0 0 ID FALSE FALSE

########
## pintar bonito las gráficas cplxCrnch_xWeighted_STAN_Summary
## Las coordenadas están en las columnas "xW_V_stan" y "xW_beta_stan" de cmplxcruncher.xlsx.
library(plotrix)

# extraemos los nombres de muestra desde df$Col1
myfrom <- 2
df <- readWorksheet(loadWorkbook("cmplxcruncher.xlsx"), sheet = 1)
rownames(df) <- unlist(lapply(strsplit(df$Col1, split="_"), FUN=function(x) paste(x[myfrom:length(x)], collapse="_")))

mycolors <- c("red", "blue", "green", "orange", "cyan")
mycols <- c(); i <- 1
mylegend <- myfill <- c()
for(g in unique(GROUPS)){
	mycols <- c(mycols, rep(mycolors[i], length(groups.l[[g]])))
	myfill <- c(myfill, mycolors[i])
	mylegend <- c(mylegend, g)
	i <- i+1
}

maxerrx <- max(df[, "xW_V_err_stan"])
maxerry <- max(df[, "xW_beta_err_stan"])

myx <- myy <- myxerr <- myyerr <- mylabels <- c()
for(g in unique(GROUPS)){
	tem <- groups.l[[g]]
	myx <- c(myx, df[tem, "xW_V_stan"])
	myy <- c(myy, df[tem, "xW_beta_stan"])
	myxerr <- c(myxerr, df[tem, "xW_V_err_stan"])
	myyerr <- c(myyerr, df[tem, "xW_beta_err_stan"])
	mylabels <- c(mylabels, tem)
}

df <- t(as.matrix(myx))
colnames(df) <- mylabels
df <- as.data.frame(df)
df$ID <- "xW_V_stan"
df <- df[, c("ID", setdiff(colnames(df), "ID"))]
write.table(df, file="xW_V_stan.tsv", sep="\t", row.names=FALSE, quote=FALSE)

df <- t(as.matrix(myy))
colnames(df) <- mylabels
df <- as.data.frame(df)
df$ID <- "xW_beta_stan"
df <- df[, c("ID", setdiff(colnames(df), "ID"))]
write.table(df, file="xW_beta_stan.tsv", sep="\t", row.names=FALSE, quote=FALSE)

df <- rbind(myx, myy)
colnames(df) <- mylabels
df <- as.data.frame(df)
df$ID <- c("xW_V_stan", "xW_beta_stan")
df <- df[, c("ID", setdiff(colnames(df), "ID"))]
write.table(df, file="xW_V_beta_stan.tsv", sep="\t", row.names=FALSE, quote=FALSE)

pdf(file="cplxCrnch_xWeighted_STAN_Summary.bonito.pdf")
plot(x=myx, y=myy, type="n", xlab="xW_V_stan", ylab="xW_beta_stan", main="Overall xWeighted_STAN cmplxcruncher Fit Summary",
	xlim=c(min(-2.01, min(myx)-maxerrx), max(2.01, max(myx)+maxerrx)), ylim=c(min(-2.01, min(myy)-maxerry), max(2.01, max(myy)+maxerry)))
draw.circle(x=0, y=0, radius=2, col=colors()[439], border="white")
draw.circle(x=0, y=0, radius=1, col=colors()[542], border=colors()[439])
segments(x0=myx-myxerr, y0=myy, x1=myx+myxerr, y1=myy, lty=3, lwd=0.7, col=colors()[210])
segments(x0=myx, y0=myy-myyerr, x1=myx, y1=myy+myyerr, lty=3, lwd=0.7, col=colors()[210])
text(x=myx, y=myy, labels=mylabels, col=mycols, cex=0.5)
legend("topleft", legend=mylegend, fill=myfill)
dev.off()

## La elección del grupo control sólo conlleva un cambio de escala que afecta a las posiciones de los puntos respecto a los círculos:
## y en realidad la distancia relativa entre puntos no cambia al cambiar de grupo control. El cambio de escala lo que hace es intentar que se queden
## el máximo número de muestras del grupo control dentro del círculo de confianza.
## Si la primera vez que se genera la gráfica V-beta weighted_STAN, aparece un grupo más estable que el que hemos designado como control, deberíamos
## ejecutar de nuevo complexCruncher redefiniendo como grupo control el que es el más estable.


