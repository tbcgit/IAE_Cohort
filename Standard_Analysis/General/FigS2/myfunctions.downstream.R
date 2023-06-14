suppressPackageStartupMessages(library(labdsv))
suppressPackageStartupMessages(library(pvclust))
suppressPackageStartupMessages(library(beeswarm))
suppressPackageStartupMessages(library(vegan))

# -------------------------------------------------------------
#
# MATRIX:
#   Abundances matrix where the columns represent samples (observations)
#   and the rows represent features to study (taxons, genes, ...).
#
# GROUPS:
#   Group assigned to each sample.
#   Named vector where the values represent distinct groups/conditions and 
#   the names represent the samples.
#
# COLOURS:
#   Colour assigned to each sample.
#   Named vector where the values represent distinct colours and 
#   the names represent the samples.
#
# LEGENDS:
#   Legend/Label assigned to each colour.  
#   Named vector where the values represent distinct labels and 
#   the names represent the colours.
#
#
# Colours are codified by name ("red", "blue", "green", ...).
# Usually colours matches to groups.
#
# -------------------------------------------------------------

doEnrich <- function(mytb, ind, ind.rest, cond, categ.descr){
	sel.go <- unlist(strsplit(mytb$categs.id[ind], split="[|]")); sel.go <- sel.go[which(!is.na(sel.go))]
	if(length(sel.go)>0){
	rest.go <- unlist(strsplit(mytb$categs.id[ind.rest], split="[|]")); rest.go <- rest.go[which(!is.na(rest.go))]
	k <- length(sel.go)
	tot <- k + length(rest.go)
	p.vals <- c(); number.in.sel <- c(); total.in.sel <- c(); number.in.rest <- c(); total.in.rest <- c()
	sel.go.un <- unique(sel.go)
	for(go in sel.go.un){
		q <- length(which(sel.go==go)); m <- q + length(which(rest.go==go))
		n <- length(which(sel.go!=go)) + length(which(rest.go!=go))
		p.vals <- c(p.vals, phyper(q=q-1, m=m, n=n, k=k, lower.tail=FALSE))
		number.in.sel <- c(number.in.sel, q)
		total.in.sel <- c(total.in.sel, k)
		number.in.rest <- c(number.in.rest, m-q)
		total.in.rest <- c(total.in.rest, tot)
	}
	ind <- which(is.element(sel.go.un, names(categ.descr)))
	res <- data.frame(categ.id=sel.go.un[ind], descr=categ.descr[sel.go.un[ind]], p.vals=p.vals[ind], number.in.sel=number.in.sel[ind],
			total.in.sel=total.in.sel[ind], number.in.rest=number.in.rest[ind], total.in.rest=total.in.rest[ind])
	res <- res[order(res$p.vals, decreasing=FALSE), ]
	write.table(res, file=paste("enrich", cond, "csv", sep="."), quote=FALSE, sep="\t", row.names=FALSE)
	}
}

diversitat <- function(taula){
    	shannon <- renyi(taula, scale=1)	# si scale=alpha=1 renyi==shannon	
	riquesa <- t(apply(taula, 1, estimateR))
	kk <- matrix(NA, dim(taula)[1], 6)
	dimnames(kk) <- list(dimnames(taula)[[1]], c("N", "Shannon", "Chao1", "SE.Chao1", "ACE", "SE.ACE"))
	kk[,1] <- riquesa[,1]
	kk[,2] <- shannon
	kk[,3:6] <- riquesa[,2:5]
	return(kk)
}

##
getCors <- function(y, M, y.name, tag, doFiles=TRUE){
	mydonors <- intersect(names(y), colnames(M))
	y <- y[mydonors]; M <- M[, mydonors]
	res.l <- vector("list", nrow(M)); names(res.l) <- rownames(M)
	for(i in 1:length(res.l)){
		dat.lm <- data.frame(y=y, x=M[i, ])
		res.lm <- lm(y~x, data=dat.lm)
		res.cor <- cor.test(x=M[i, ], y=y, method="spearman")
		res.l[[i]] <- list(res.lm=res.lm, dat.lm=dat.lm, cor=res.cor$estimate, cor.p.val=res.cor$p.value)
	}
	cor.v <- unlist(lapply(res.l, FUN=function(x) x$cor))
	cor.p.val <- unlist(lapply(res.l, FUN=function(x) x$cor.p.val))
	lm.p.val <- unlist(lapply(res.l, FUN=function(x) summary(x$res.lm)$coefficients[2,4]))
	mytb <- data.frame(y=rep(y.name, length(res.l)), x=names(res.l), spearman.cor=cor.v,
			cor.p.val=cor.p.val, lm.p.val=lm.p.val, x.descr=DESCRIPTION[names(res.l)], stringsAsFactors=FALSE)
	mytb <- mytb[order(abs(mytb$spearman.cor), decreasing=TRUE), ]
	mytb$cor.adj.p.val <- p.adjust(mytb$cor.p.val, method="fdr")
	if(doFiles){
	write.table(mytb, file=paste("corrs", tag, "tsv", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
	cor.p.val.cut <- 0.1
	ind <- which(mytb$cor.p.val<=cor.p.val.cut)
	if(length(ind)>1){
		mytb.cut <- mytb[ind, ]
		res.l <- res.l[mytb.cut$x]
		pdf(file=paste("corrs", tag, "pdf", sep="."))
		for(i in 1:nrow(mytb.cut)){
			x <- res.l[[i]]$dat.lm$x; y <- res.l[[i]]$dat.lm$y
			plot(x=x, y=y, type="n", xlab=mytb.cut[i, "x"], ylab=mytb.cut[i, "y"], cex.lab=0.6)
			text(x=x, y=y, labels=mydonors, col="blue", cex=0.6)
			abline(res.l[[i]]$res.lm, col="red")
			title(sub=paste("cor:", format(mytb.cut[i, "spearman.cor"], digits=2),
						"cor.p.val:", format(mytb.cut[i, "cor.p.val"], digits=2),
						"lm.p.val:", format(mytb.cut[i, "lm.p.val"], digits=2), sep=" "),
				main=DESCRIPTION[mytb.cut[i, "x"]], cex.main=0.4, cex.sub=0.6)
		}
		dev.off()
	}
	}
	return(mytb)
}

getCorsResponseA <- function(y, M, y.name, tag, doFiles=TRUE){
	mydonors <- intersect(names(y), colnames(M))
	y <- y[mydonors]; M <- M[, mydonors]
	res.l <- vector("list", nrow(M)); names(res.l) <- rownames(M)
	for(i in 1:length(res.l)){
		dat.lm <- data.frame(y=y, x=M[i, ])
		res.lm <- lm(y~x, data=dat.lm)
		res.cor <- cor.test(x=M[i, ], y=y, method="spearman")
		res.l[[i]] <- list(res.lm=res.lm, dat.lm=dat.lm, cor=res.cor$estimate, cor.p.val=res.cor$p.value)
	}
	cor.v <- unlist(lapply(res.l, FUN=function(x) x$cor))
	cor.p.val <- unlist(lapply(res.l, FUN=function(x) x$cor.p.val))
	lm.p.val <- unlist(lapply(res.l, FUN=function(x) summary(x$res.lm)$coefficients[2,4]))
	mytb <- data.frame(y=rep(y.name, length(res.l)), x=names(res.l), spearman.cor=cor.v,
			cor.p.val=cor.p.val, lm.p.val=lm.p.val, x.descr=DESCRIPTION[names(res.l)], stringsAsFactors=FALSE)
	mytb <- mytb[order(abs(mytb$spearman.cor), decreasing=TRUE), ]
	mytb$cor.adj.p.val <- p.adjust(mytb$cor.p.val, method="fdr")
	if(doFiles){
	write.table(mytb, file=paste("corrs", tag, "tsv", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
	cor.p.val.cut <- 0.1
	ind <- which(mytb$cor.p.val<=cor.p.val.cut)
	if(length(ind)>1){
		mytb.cut <- mytb[ind, ]
		res.l <- res.l[mytb.cut$x]
		pdf(file=paste("corrs", tag, "pdf", sep="."))
		for(i in 1:nrow(mytb.cut)){
			x <- res.l[[i]]$dat.lm$x; y <- res.l[[i]]$dat.lm$y
			plot(x=x, y=y, type="n", xlab=mytb.cut[i, "x"], ylab=mytb.cut[i, "y"], cex.lab=0.6)
			text(x=x, y=y, labels=mydonors, col="blue", cex=0.6)
			abline(res.l[[i]]$res.lm, col="red")
			title(sub=paste("cor:", format(mytb.cut[i, "spearman.cor"], digits=2),
						"cor.p.val:", format(mytb.cut[i, "cor.p.val"], digits=2),
						"lm.p.val:", format(mytb.cut[i, "lm.p.val"], digits=2), sep=" "),
				main=DESCRIPTION[mytb.cut[i, "x"]], cex.main=0.4, cex.sub=0.6)
		}
		dev.off()
	}
	}
	return(mytb)
}

getCorsResponseB <- function(y, M, y.name, tag, doFiles=TRUE){
	mydonors <- intersect(names(y), colnames(M))
	y <- y[mydonors]; M <- M[, mydonors]
	res.l <- vector("list", nrow(M)); names(res.l) <- rownames(M)
	for(i in 1:length(res.l)){
		dat.lm <- data.frame(y=M[i, ], x=y)
		res.lm <- lm(y~x, data=dat.lm)
		res.cor <- cor.test(x=y, y=M[i, ], method="spearman")
		res.l[[i]] <- list(res.lm=res.lm, dat.lm=dat.lm, cor=res.cor$estimate, cor.p.val=res.cor$p.value)
	}
	cor.v <- unlist(lapply(res.l, FUN=function(x) x$cor))
	cor.p.val <- unlist(lapply(res.l, FUN=function(x) x$cor.p.val))
	lm.p.val <- unlist(lapply(res.l, FUN=function(x) summary(x$res.lm)$coefficients[2,4]))
	mytb <- data.frame(y=names(res.l), x=rep(y.name, length(res.l)), spearman.cor=cor.v,
			cor.p.val=cor.p.val, lm.p.val=lm.p.val, y.descr=DESCRIPTION[names(res.l)], stringsAsFactors=FALSE)
	mytb <- mytb[order(abs(mytb$spearman.cor), decreasing=TRUE), ]
	mytb$cor.adj.p.val <- p.adjust(mytb$cor.p.val, method="fdr")
	if(doFiles){
	write.table(mytb, file=paste("corrs", tag, "tsv", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
	cor.p.val.cut <- 0.1
	ind <- which(mytb$cor.p.val<=cor.p.val.cut)
	if(length(ind)>1){
		mytb.cut <- mytb[ind, ]
		res.l <- res.l[mytb.cut$y]
		pdf(file=paste("corrs", tag, "pdf", sep="."))
		for(i in 1:nrow(mytb.cut)){
			x <- res.l[[i]]$dat.lm$x; y <- res.l[[i]]$dat.lm$y
			plot(x=x, y=y, type="n", xlab=mytb.cut[i, "x"], ylab=mytb.cut[i, "y"], cex.lab=0.6)
			text(x=x, y=y, labels=mydonors, col="blue", cex=0.6)
			abline(res.l[[i]]$res.lm, col="red")
			title(sub=paste("cor:", format(mytb.cut[i, "spearman.cor"], digits=2),
						"cor.p.val:", format(mytb.cut[i, "cor.p.val"], digits=2),
						"lm.p.val:", format(mytb.cut[i, "lm.p.val"], digits=2), sep=" "),
				main=DESCRIPTION[mytb.cut[i, "y"]], cex.main=0.4, cex.sub=0.6)
		}
		dev.off()
	}
	}
	return(mytb)
}

doFiltering <- function(dat, cutoff, group, mypct){
	# dat ha de tener muestras en columnas y variables en rownames
	good <- c()
	for(mygroup in unique(as.character(group))){
		ind <- which(group==mygroup); n <- round(length(ind)*mypct)
		if(n<3) n <- 3
		if(n>length(ind)) n <- length(ind)
		dat.temp <- dat[, ind]
		cond <- apply(dat.temp, MARGIN=1, FUN=function(x) if(length(which(abs(x)>cutoff))>=n) return(TRUE) else return(FALSE))
		good <- c(good, rownames(dat.temp)[which(cond==TRUE)])
	}
	return(unique(good))
}

doFiltering2 <- function(dat, cutoff, group, mynum=5){
	# dat ha de tener muestras en columnas y variables en rownames
	good <- c()
	for(mygroup in unique(as.character(group))){
		dat.temp <- dat[, which(group==mygroup)]
		cond <- apply(dat.temp, MARGIN=1, FUN=function(x) if(length(which(x!=cutoff))>=mynum) return(TRUE) else return(FALSE))
		good <- c(good, rownames(dat.temp)[which(cond==TRUE)])
	}
	return(unique(good))
}

################################################
# Remove a feature (row) in the matrix 
# if the value of a specific number of observations (N) in each group
# is not greater than a specific value (CUTOFF)
################################################
remove.lowvalues.bygroup <- function(MATRIX, GROUPS, CUTOFF=0, N=3)
{
  # Distinct groups or conditions
  groups.unique <- unique(GROUPS)

  rows.ok <- sapply(groups.unique, function(group)
  {
    # Observations associated to the group
    cols.group <- names(GROUPS)[which(GROUPS==group)]
    matrix.group <- MATRIX[,cols.group]
    
    # Valid rows in the group
    rows.group.ok <- apply(matrix.group, MARGIN=1, function(row.values.group)
    {
      return(length(which(row.values.group > CUTOFF)) >= N)
    })
    
    return(which(rows.group.ok==TRUE))
  })
  
  # At least the row is valid in one group
  rows.ok <- Reduce(union, rows.ok)
  
  MATRIX <- MATRIX[rows.ok,]

  return(MATRIX)
}


#################################################
# Hierarchical clustering of the observations
#################################################
hierarchical.clustering <- function(MATRIX, FILE="")
{
  # Hierarchical clustering via multiscale bootstrap resampling
  matrix.pvclust <- pvclust(MATRIX, method.hclust="complete", method.dist="euclidean", nboot=200, quiet=TRUE)
  
  if(FILE!="") pdf(FILE)

  # Dendogram with p-values
  plot(matrix.pvclust, main="")
  
  # Add rectangles around groups highly supported by the data
  pvrect(matrix.pvclust, alpha=0.95)
  
  if(FILE!="") temp <- dev.off()
}


#################################################
# Test differences between groups of observations (columns)
# for each feature (row) in the matrix.
#
# TEST: {t.test, wilcox.test}
#################################################
differential.tests <- function(MATRIX, GROUPS, GROUP1, GROUP2, TEST, PAIRED, CUTOFF, COLNAME, FILE="", DESCRIPTION=NULL,
						myfitType="parametric", meta=NULL, datnorm=NULL)
{

	if(length(apropos("mynorm"))==0) mynorm <- TRUE
  # Observations in each group
  items.group1 <- names(GROUPS)[GROUPS==GROUP1]
  items.group2 <- names(GROUPS)[GROUPS==GROUP2]
  
  # Differences between both groups
  if(TEST=="t.test"|TEST=="wilcox.test"){
  p.values <- apply(MATRIX, MARGIN=1, function(values) 
  {
    do.call(TEST, list(x=values[items.group1], y=values[items.group2], paired=PAIRED))$p.value
  })
  }
  
	#if(TEST=="DESeq.test"){
	#	countTable <- MATRIX
	#	if(any(countTable==0)){
			#countTable <- countTable + matrix(round(runif(ncol(countTable)*nrow(countTable), 1, 10)), ncol=ncol(countTable), nrow(countTable))
	#		countTable <- countTable + 1
	#	}
	#	condition <- as.factor(GROUPS)
	#	cds <- newCountDataSet(countTable, condition)
	#	cds <- estimateSizeFactors(cds)
	#	if(PAIRED==FALSE){
	#		cds <- estimateDispersions(cds, fitType=myfitType)
	#		res <- nbinomTest(cds, GROUP1, GROUP2)
	#		p.values <- res$pval
	#	} else {
	#		cds <- estimateDispersions(cds, method="blind", fitType=myfitType)
	#		vsd <- getVarianceStabilizedData(cds)
	#		if(mynorm==TRUE) vsd <- prop.table(vsd, margin=2)*100
	#		p.values <- apply(vsd, MARGIN=1, function(values) 
  	#		{
    	#			do.call("wilcox.test", list(x=values[items.group1], y=values[items.group2], paired=PAIRED))$p.value
  	#		})
	#	}
	#}

	if(TEST=="DESeq2.test"){
		if(PAIRED==FALSE){
			coldata <- data.frame(sample=names(GROUPS), condition=as.factor(GROUPS))
			rownames(coldata) <- names(GROUPS)
			dds <- DESeqDataSetFromMatrix(countData = MATRIX, colData = coldata, design = ~ condition)
		} else {
			coldata <- data.frame(sample=names(GROUPS), condition=as.factor(GROUPS), donor=as.factor(meta$donor))
			rownames(coldata) <- names(GROUPS)
			dds <- DESeqDataSetFromMatrix(countData = MATRIX, colData = coldata, design = ~ donor + condition)
		}
		dds <- DESeq(dds, fitType=myfitType)
		res <- lfcShrink(dds, coef=2, type="apeglm")
		res <- as.data.frame(res@listData)
		p.values <- res$pvalue
		p.values.adj <- res$padj
		#FC <- res$log2FoldChange
	}

  if(mynorm==TRUE&!is.null(datnorm)&TEST=="DESeq2.test") MATRIX <- datnorm
  if(TEST!="DESeq2.test"){
	p.values.adj <- p.adjust(p.values, method="fdr")
  }
  
  # Mean value for each group
  means.group1 <- apply(MATRIX, MARGIN=1, function(values) {mean(values[items.group1])})
  means.group2 <- apply(MATRIX, MARGIN=1, function(values) {mean(values[items.group2])})

	#if(TEST!="DESeq2.test"){
		mymin1 <- min(means.group1[which(means.group1!=0)])
		mymin2 <- min(means.group2[which(means.group2!=0)])
		mymin <- mean(c(mymin1, mymin2))
		FC <- log2((means.group1+mymin)/(means.group2+mymin))
	#}

  df <- data.frame(p.value=p.values, adj.p.value=p.values.adj, 
                   mean1=means.group1, mean2=means.group2, log2FC=FC, stringsAsFactors=FALSE)

  # Lower p.value --> higher difference
  df <- df[order(df$p.value, decreasing=FALSE),]
  
  # Get only significant differences
  df2plot <- df[which(df$adj.p.value < CUTOFF),]
  if(nrow(df2plot)<5) df2plot <- df[which(df$p.value < CUTOFF),]
  
  # Write p.values into a text file
  df.aux <- cbind(rownames(df), df)
  colnames(df.aux)[1] <- COLNAME
  if(!is.null(DESCRIPTION)) df.aux$DESCRIPTION <- DESCRIPTION[as.character(df.aux[, 1])]

  if(FILE!="") write.table(df.aux, paste(FILE,"tsv",sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  # Boxplot values of rows with significant differences
  if(FILE!=""){

  pdf(file=paste(FILE,"pdf",sep="."))
  for(sign.row in rownames(df2plot))
  {
    values <- list(MATRIX[sign.row,items.group1], MATRIX[sign.row,items.group2])
    names(values) <- c(GROUP1, GROUP2)
    if(!is.null(DESCRIPTION)) mymain <- DESCRIPTION[sign.row] else mymain <- ""
    beeswarm(values, pch=16, col=c("red","blue"), las=1, xlab="Condition", ylab=sign.row, main=mymain, cex.lab=0.6, cex.main=0.5)
    boxplot (values, las=1,  add=TRUE, col=NULL)
    legend("topright", legend=c(paste("p.value:"         , format(df2plot[sign.row,"p.value"    ],digits=2)),
                                paste("adjusted p.value:", format(df2plot[sign.row,"adj.p.value"],digits=2))), cex=0.7)
  }
  temp <- dev.off()
  
  }

  return(df.aux)
}


#################################################
# Assign the same colour to each observation in the same group
# and a distinct colour for distinct groups.
#################################################
get.colours <- function(GROUPS)
{
  mypallete <- c("red", "blue", "green3", "gold4", "darkviolet",
                 "orange", "cyan3", "brown", "cornflowerblue", "coral", 
                 "turquoise", "slategray", "yellowgreen", "plum")
  
  groups.unique <- unique(GROUPS)
  groups.factor <- factor(GROUPS, ordered=TRUE, levels=groups.unique)

  groups.numeric <- as.numeric(groups.factor)
  
  if(length(groups.unique) > length(mypallete))
    stop("Too many groups, add colours to the pallete.")
  
  mycolours <- mypallete[groups.numeric]
  names(mycolours) <- names(GROUPS)
  
  mylegends        <- groups.unique
  names(mylegends) <- mypallete[unique(groups.numeric)]

  return(list(colour=mycolours, legend=mylegends))
}


#################################################
# Principal coordinates analysis of distances between observations
#
# Moves each observation (column)
# to a new (eigen vectors) coordinates system.
# Possible correlations between features (rows) are lost
# in the new coordinates system.
#################################################
principal.coordinates.analysis <- function(MATRIX, GROUPS, TITLE="", METHOD="euclidean", ELLIPSE=FALSE, FILE="", plot.arrows=FALSE, sum2lims=0.15)
{
  if(class(MATRIX)=="dist"){
	matrix.dist <- MATRIX
  } else {
  	matrix.dist <- c()
  	# Distances/Dissimilarities between row vectors
  	# matrix.dist <- labdsv::dsvdis(t(MATRIX), METHOD)
  	matrix.dist <- vegdist(t(MATRIX), METHOD)
  }

  # Principal coordinates analysis
  # 3 first principal components
  matrix.dist.pco <- pco(matrix.dist, k=3)

  if(FILE!="") pdf(FILE)

  # The greatest variance lies on the first coordinate, the second highest variance lies on the second coordinate and so on.
  # The proportion of the variance that each eigenvector represents can be calculated by dividing the eigenvalue corresponding to that eigenvector by the sum of all eigenvalues
  var1 <- 100*matrix.dist.pco$eig[1]/sum(matrix.dist.pco$eig)
  var2 <- 100*matrix.dist.pco$eig[2]/sum(matrix.dist.pco$eig)
  var3 <- 100*matrix.dist.pco$eig[3]/sum(matrix.dist.pco$eig)

  myxlab <- paste("PCO1: ", format(var1, digits=4), "%", sep="")
  myylab <- paste("PCO2: ", format(var2, digits=4), "%", sep="")
  myzlab <- paste("PCO3: ", format(var3, digits=4), "%", sep="")

  # Coordinates of each observation in the new system
  vx <- matrix.dist.pco$points[,1]
  vy <- matrix.dist.pco$points[,2]
  vz <- matrix.dist.pco$points[,3]
  
  # almacenamos coordenadas de los puntos que aparecerán en el plot
	write.csv(matrix.dist.pco$points, file=paste(sub(".pdf", "", FILE), "csv", sep="."))

  # Colour for each observation
  #GROUPS <- GROUPS[colnames(MATRIX)]
  mycolours <- get.colours(GROUPS)

  # Plot observations
  layout(matrix(data=c(1,2), nrow=1, ncol=2), 
         widths=c(10,3), heights=c(1,1))

  xlim <- c(min(vx)-sum2lims, max(vx)+sum2lims)
  ylim <- c(min(vy)-sum2lims, max(vy)+sum2lims)

  par(mar=c(5,5,2,0)) # Bottom, Left, Top, Right
  plot(x=vx, y=vy, xlab=myxlab, ylab=myylab, las=1, main=TITLE, type="n", xlim=xlim, ylim=ylim)
  #text(x=vx, y=vy, labels=names(GROUPS), col=mycolours$colour, cex=0.5)
  points(x=vx, y=vy, pch=20, col=mycolours$colour)

  if(plot.arrows==TRUE){
	arrows(x0=rep(0, length(vx)), y0=rep(0, length(vy)), x1=vx, y1=vy, col=mycolours$colour, length = 0.0)
  }

  if(ELLIPSE)
  {
    matrix.dist.pco.ellipse <- ordiellipse(matrix.dist.pco,
                draw="polygon", alpha=0.2, kind="sd",
                groups=factor(GROUPS, ordered=TRUE, levels=mycolours$legend),
                col=names(mycolours$legend), border=names(mycolours$legend))
  }
  
  par(mar=c(5,1,1.5,0))
  plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  legend("topleft", mycolours$legend, pch=15, col=names(mycolours$legend), bty="n", cex=0.7)

  if(FILE!="") temp <- dev.off()
}


#################################################
# Test differences between the groups of observations
# for groups of related features (multivariate analysis of variance).
# Extracts major gradients among combinations of explanatory features.
#################################################
correspondence.analysis <- function(MATRIX, GROUPS, COLOURS=c(), LEGENDS=c(), FILE="",
						plot.centroids=TRUE, plot.sites.labels=TRUE, plot.arrows=TRUE)
{
  # Distinct groups or conditions
  groups.unique <- unique(GROUPS)
  groups.unique.aux <- paste("group", groups.unique, sep="")
  names(groups.unique.aux) <- groups.unique
  groups.factor <- factor(GROUPS, ordered=TRUE, levels=groups.unique)

  # Data frame of observations with groups
  matrix.df <- as.data.frame(t(MATRIX), stringAsFactors=FALSE)
  matrix.groups.df <- cbind(matrix.df, group=as.factor(GROUPS))
  
  # Canonical correspondence analysis
  matrix.cca <- vegan::cca(matrix.df ~ group, data=matrix.groups.df)
  matrix.cca.summary <- summary(matrix.cca)
  matrix.cca.anova <- anova(matrix.cca)
  matrix.cca.p.value <- format(matrix.cca.anova[["Pr(>F)"]][1], digits=2)
  
  # Analysis of variance using distance matrices
  matrix.adonis <- adonis(matrix.df ~ group, data=matrix.groups.df)
  matrix.adonis.p.value <- format(matrix.adonis$aov.tab["group", "Pr(>F)"], digits=2)

  
  if(length(groups.unique) > 2) {p.var <- matrix.cca.summary$concont$importance["Proportion Explained",c("CCA1","CCA2")]
  } else                        {p.var <- matrix.cca.summary$cont   $importance["Proportion Explained",c("CCA1","CA1" )]}
  
  p.var <- format(p.var*100, digits=4)
  component1 <- names(p.var)[1]
  component2 <- names(p.var)[2]
  
  myxlab <- paste(component1, paste(p.var[component1],"%", sep=""))
  myylab <- paste(component2, paste(p.var[component2],"%", sep=""))

  if(FILE!="") pdf(FILE)

  mytitle <- paste("CCA p-value: ", matrix.cca.p.value, " - ADONIS p-value: ", matrix.adonis.p.value, sep="")
  matrix.cca.plot <- plot(matrix.cca, las=1, type="none", xlab=myxlab, ylab=myylab, font.main=1, main=mytitle)

	# almacenamos coordenadas de los puntos que aparecerán en el plot
	write.csv(rbind(matrix.cca.plot$centroids, matrix.cca.plot$sites), file=paste(sub(".pdf", "", FILE), "csv", sep="."))

  if(length(COLOURS)==0)
  {
    mycolours <- get.colours(GROUPS)
    COLOURS <- mycolours$colour
  }

  mycolour.groups <- c()
  
  for(mygroup in groups.unique)
  {
    items.group <- names(GROUPS)[GROUPS==mygroup]
    
    # Colour for observations in the group
    mycolour.group  <- unique(COLOURS[items.group])
    mycolour.groups <- c(mycolour.groups, mycolour.group)
    names(mycolour.groups)[length(mycolour.groups)] <- mygroup
    
    mygroup.aux <- groups.unique.aux[mygroup]

    # Line between the origin and the center of the ellipse (group)
    if(plot.arrows==TRUE){
    arrows(x0=0, y0=0, 
           x1=matrix.cca.plot$centroids[mygroup.aux, component1], 
           y1=matrix.cca.plot$centroids[mygroup.aux, component2], 
           col=mycolour.group, length=0)
    }

    # Name of the group in the center of the ellipse
    if(plot.centroids==TRUE){
    text(labels=mygroup,
         x=matrix.cca.plot$centroids[mygroup.aux, component1], 
         y=matrix.cca.plot$centroids[mygroup.aux, component2],
         col=mycolour.group, cex=0.7)
    }

    # Observations in the group
    if(plot.sites.labels==TRUE){
    text(x=matrix.cca.plot$sites[items.group,component1], 
         y=matrix.cca.plot$sites[items.group,component2], 
         col=mycolour.group, cex=0.4, labels=items.group)
    } else {
	points(x=matrix.cca.plot$sites[items.group,component1], 
         y=matrix.cca.plot$sites[items.group,component2], 
         col=mycolour.group, pch=20)
    }
  }
  
  # Ellipse for each group
  ordiellipse(matrix.cca, groups=groups.factor, col=mycolour.groups, alpha=0.2, draw="polygon", kind="sd", border=mycolour.groups)
  
  if(length(LEGENDS) > 0)
    legend("topright", LEGENDS, pch=15, col=names(LEGENDS))
  
  if(FILE!="") temp <- dev.off()
}

doBeeswarm <- function(MATRIX, GROUPS, FILE="", p.value=NULL, adj.p.value=NULL, test="", CUTOFF=0.05){
	colors <- c("red", "blue", "orange", "cyan")
	if(test=="kruskal.test"){
		p.value <- c()
		for(i in 1:nrow(MATRIX)){
			dat.temp <- data.frame(y=MATRIX[i, ], group=GROUPS)
			p.value <- c(p.value, kruskal.test(y ~ group, data=dat.temp)$p.value)
		}
		names(p.value) <- rownames(MATRIX)
		adj.p.value <- p.adjust(p.value, method="fdr")
		mytb <- data.frame(taxa=names(p.value), p.value=p.value, adj.p.value=adj.p.value)
		myord <- order(mytb$p.value, decreasing=FALSE)
		mytb <- mytb[myord, ]
		p.value <- p.value[myord]; adj.p.value <- adj.p.value[myord]
		ind <- which(mytb$p.value<=CUTOFF)
		write.table(mytb, paste(FILE,"tsv",sep="."), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
		vars2plot <- as.character(mytb$taxa)[ind]
		MATRIX <- MATRIX[vars2plot, ]
	}
	pdf(file=paste(FILE,"pdf",sep="."))
	for(i in 1:nrow(MATRIX)){
		dat.temp <- data.frame(y=MATRIX[i, ], group=GROUPS)
		beeswarm(y ~ group, data=dat.temp, pch=16, col=colors, las=1, xlab="Condition", ylab=rownames(MATRIX)[i], main="", cex.lab=0.6)
		boxplot (y ~ group, data=dat.temp, las=1, add=TRUE, col=NULL)
		legend("topright", legend=c(paste("p.value:"         , format(p.value[i], digits=2)),
                                paste("adjusted p.value:", format(adj.p.value[i], digits=2))), cex=0.7)
	}
	dev.off()
	return(mytb)
}

getVars <- function(dat, GROUPS, group, cut.up=0, cut.down=0, n, direc){
	ind <- which(GROUPS==group)
	if(n>length(ind)) n <- length(ind)
	if(direc=="up") cond <- apply(dat[, ind], MARGIN=1, FUN=function(x) if(length(which(x>cut.up))>=n&mean(x)>cut.up) return(TRUE) else return(FALSE))
	if(direc=="down") cond <- apply(dat[, ind], MARGIN=1, FUN=function(x) if(length(which(x<cut.down))>=n&mean(x)<cut.down) return(TRUE) else return(FALSE))
	return(rownames(dat)[which(cond==TRUE)])
}

doVenn3 <- function(GROUP1, GROUP2, GROUP3, GROUP1.vals, GROUP2.vals, GROUP3.vals, tag){
	mycomp <- paste(GROUP1, "vs", GROUP2, "vs", GROUP3, sep=".") 
	groups <- list(GROUP1.vals, GROUP2.vals, GROUP3.vals)
	names(groups) <- c(GROUP1, GROUP2, GROUP3)
	inter <- intersect(intersect(groups[[1]], groups[[2]]), groups[[3]])
	inter.1.2 <- setdiff(intersect(groups[[1]], groups[[2]]), inter)
	inter.1.3 <- setdiff(intersect(groups[[1]], groups[[3]]), inter)
	inter.2.3 <- setdiff(intersect(groups[[2]], groups[[3]]), inter)
	only.1 <- setdiff(groups[[1]], unique(c(groups[[2]], groups[[3]])))
	only.2 <- setdiff(groups[[2]], unique(c(groups[[1]], groups[[3]])))
	only.3 <- setdiff(groups[[3]], unique(c(groups[[1]], groups[[2]])))
	write(inter, file=paste("intersection", tag, "txt", sep="."), sep="\n")
	write(inter.1.2, file=paste("intersection", GROUP1, GROUP2, "txt", sep="."), sep="\n")
	write(inter.1.3, file=paste("intersection", GROUP1, GROUP3, "txt", sep="."), sep="\n")
	write(inter.2.3, file=paste("intersection", GROUP2, GROUP3, "txt", sep="."), sep="\n")
	write(only.1, file=paste("only.group", GROUP1, "txt", sep="."), sep="\n")
	write(only.2, file=paste("only.group", GROUP2, "txt", sep="."), sep="\n")
	write(only.3, file=paste("only.group", GROUP3, "txt", sep="."), sep="\n")
	pdf(file=paste("venn.diagram", mycomp, "pdf", sep="."))
	venn(groups)
	dev.off()
}

plotPaired <- function(dat, feature, group, group.x, group.y, individual.factor, DESCRIPTION, p.val=NULL){
	mycols <- c("blue", "red")
	ind.x <- which(group==group.x)
	ind.y <- which(group==group.y)
	individual.factor.v <- as.character(individual.factor)
	names(individual.factor.v) <- colnames(dat)
	x.v <- as.vector(as.matrix(dat[feature, ind.x])); y.v <- as.vector(as.matrix(dat[feature, ind.y]))
	x.v <- x.v[order(names(individual.factor.v)[ind.x])]
	y.v <- y.v[order(names(individual.factor.v)[ind.y])]
	names(x.v) <- names(y.v) <- individual.factor.v[order(names(individual.factor.v)[ind.x])]
	if(is.null(p.val)) p.val <- wilcox.test(x=x.v, y=y.v, paired=TRUE)$p.value
	v <- c(x.v, y.v)
	group <- as.character(group)[c(ind.x, ind.y)]
	plot(x=1:length(v), y=v, type="n", xaxt="n", ylab="", xlab="", xlim=c(0,8))
	points(x=rep(2, length(x.v)), y=x.v, col=mycols[1])
	points(x=rep(6, length(y.v)), y=y.v, col=mycols[2])
	text(x=rep(2, length(x.v)), y=x.v, labels=names(x.v), cex=0.6, col=mycols[1], pos=2)
	text(x=rep(6, length(y.v)), y=y.v, labels=names(y.v), cex=0.6, col=mycols[2], pos=4)
	axis(side=1, at=c(2,6), labels=c(group.x, group.y))
	segments(x0=rep(2, length(x.v)), y0=x.v, x1=rep(6, length(y.v)), y1=y.v, lty=2, col="grey")
	values <- list(x.v, y.v)
	names(values) <- c(group.x, group.y)
	graphics::boxplot(values, add = TRUE, at=c(2,6), col=NULL)
	title(sub=paste("p-value:", format(p.val, digits=2), sep=" "), cex.sub=0.6, cex.lab=0.6,
		ylab=feature, main=DESCRIPTION[feature], cex.main=0.5)
}

doSummary <- function(file=NULL, mytb=NULL, mincutoff=0.6, sampleColumn="cluster", groups2plot=NULL,
						GROUPS=NULL, tag=NULL, mybottom=5, samplesORD=NULL, fixFileName=NULL, ind.taxa="all", mytaxlevel="species"){
	if(is.null(mytb)){
		mytb <- read.table(file=file, sep="\t", quote="", stringsAsFactors=FALSE, header=TRUE, comment.char="")
		file <- sub("[.]tsv", "", file)
	}
	if(is.null(file)) file <- ""
	if(ind.taxa=="all") ind.taxa <- 1:nrow(mytb) else ind.taxa <- grep("Bacteria", mytb$feature)
	ind.func <- setdiff(1:nrow(mytb), ind.taxa)
	if(length(ind.taxa)>0){
		taxlevel <- mytaxlevel
		all.data <- mytb[ind.taxa, c("feature", sampleColumn)]
		all.data$sample <- as.character(all.data[, sampleColumn])
		if(!is.null(groups2plot)&!is.null(GROUPS)){
			mysamples <- c()
			for(i in 1:length(groups2plot)) mysamples <- c(mysamples, names(GROUPS)[which(GROUPS==groups2plot[i])]) 
			all.data <- all.data[which(is.element(all.data$sample, mysamples)), ]
			if(!is.null(fixFileName)) file <- fixFileName else file <- paste(groups2plot, collapse=".vs.")
		} else {
			file <- tag
		}
		if(mytaxlevel=="species") all.data$species <- unlist(lapply(strsplit(all.data$feature, split=";"), FUN=function(x) paste(x[c(2, 4:6)], collapse="|")))
		if(mytaxlevel=="genus") all.data$genus <- unlist(lapply(strsplit(all.data$feature, split=";"), FUN=function(x) paste(x[c(2, 4:5)], collapse="|")))
		all.data$phylum <- unlist(lapply(strsplit(all.data$feature, split=";"), FUN=function(x) x[1]))
		data <- table(all.data$sample, all.data[, taxlevel])
		data.pct <- prop.table(data, margin = 1) * 100
		if(!is.null(groups2plot)&!is.null(GROUPS)) data.pct <- data.pct[mysamples, ]
		# buscamos los minoritarios
		cond <- apply(data.pct, MARGIN=2, FUN=function(x) if(all(x<mincutoff)) return(TRUE) else return(FALSE))
		minor <- colnames(data.pct)[which(cond==TRUE)]
		# recalculamos table para el barplot
		if(length(minor)>0){
			for(taxa in minor){
				ind <- which(all.data[, taxlevel]==taxa)
				all.data[ind, taxlevel] <- "other"
				all.data[ind, "phylum"] <- "other"
				print(taxa)
			}
		}
		ind <- which(!duplicated(all.data[, taxlevel]))
		alltaxa <- all.data[ind, taxlevel]
		allphylum <- all.data[ind, "phylum"]
		names(allphylum) <- alltaxa
		ta.dfm <- table(all.data$sample, all.data[, taxlevel])
		ta.dfm.pct <- prop.table(ta.dfm, margin = 1) * 100
		if(!is.null(groups2plot)&!is.null(GROUPS)) ta.dfm.pct <- ta.dfm.pct[mysamples, ]
		myphylum <- allphylum[colnames(ta.dfm.pct)]
		if(length(minor)>0){
			ta.dfm.pct <- cbind(ta.dfm.pct[, setdiff(colnames(ta.dfm.pct), "other")], other=ta.dfm.pct[, "other"])
			myphylum <- c(myphylum[which(myphylum!="other")], "other")
		}
		names(myphylum) <- NULL
		if(!is.null(samplesORD)){
			ta.dfm.pct <- ta.dfm.pct[samplesORD, ]
		} else {
			ta.dfm.pct <- ta.dfm.pct[order(rownames(ta.dfm.pct), decreasing=TRUE), ]
		}
		ta.dfm.pct <- ta.dfm.pct[nrow(ta.dfm.pct):1, ]
		write.csv2(t(ta.dfm.pct), file=paste("summary.taxa", file, "csv", sep="."), quote=FALSE)
		myphylum <- sub("p[.]", "", myphylum)
		save(myphylum, file=paste("phylum", file, "RData", sep="."))
		pdf(file=paste("summary.taxa", file, "pdf", sep="."))
		mybarplot(myphylum=myphylum, ta.dfm.pct=ta.dfm.pct, "", mybottom=mybottom)
		dev.off()
	}
	if(length(ind.func)>0){
		all.data <- mytb[ind.func, c("Route_dsc", sampleColumn)]
		all.data$sample <- as.character(all.data[, sampleColumn])
		data <- table(all.data$sample, all.data$Route_dsc)
		data <- data[order(rownames(data), decreasing=TRUE), ]
		data.pct <- prop.table(data, margin = 1) * 100
		# buscamos los minoritarios
		cond <- apply(data.pct, MARGIN=2, FUN=function(x) if(all(x<mincutoff)) return(TRUE) else return(FALSE))
		minor <- colnames(data.pct)[which(cond==TRUE)]
		# recalculamos table para el barplot
		if(length(minor)>0){
			for(taxa in minor){
				ind <- which(all.data$Route_dsc==taxa)
				all.data[ind, "Route_dsc"] <- "other"
				print(taxa)
			}
		}
		ta.dfm <- table(all.data$sample, all.data$Route_dsc)
		ta.dfm.pct <- prop.table(ta.dfm, margin = 1) * 100
		pdf(file=paste("summary.func", file, "pdf", sep="."))
		for(i in 1:nrow(ta.dfm.pct)){
			my.v <- ta.dfm.pct[i, ]
			par(mar=c(3.7, 17, 1, 2) + 0.1)
			barplot(my.v, horiz=TRUE, las=2, cex.names=0.4, col=rainbow(length(my.v)), xlab="% abundance",
				main=paste(sampleColumn, rownames(ta.dfm.pct)[i], sep=" "), xlim=c(0, max(ta.dfm.pct)+1))
		}
		dev.off()
	}
}

mybarplot <- function(myphylum, ta.dfm.pct, mymain, mybottom=5, splitLegend=FALSE){	
	
	my.extreme.cols1 <- c("palevioletred4", "snow2", "rosybrown4", colors()[560], "orange4", "snow2",
					colors()[470], colors()[468], "lightgoldenrod4", "snow2",
					colors()[544], colors()[542], colors()[256], "snow2", "lightslateblue", "snow2")

	if(is.element("other", colnames(ta.dfm.pct))) other.pct <- ta.dfm.pct[, "other"]
	my.ta.dfm.pct <- ta.dfm.pct[, setdiff(colnames(ta.dfm.pct), "other")]	
	my.phylums <- myphylum[which(myphylum!="other")]
	myord <- order(my.phylums)
	my.phylums <- my.phylums[myord]
	my.ta.dfm.pct <- my.ta.dfm.pct[, myord]
	my.taxas <- colnames(my.ta.dfm.pct)
	
	my.unique.phy <- unique(my.phylums)
	my.col.lens <- c()
	for(i in 1:length(my.unique.phy)) my.col.lens <- c(my.col.lens, length(grep(my.unique.phy[i], my.phylums)))
	my.cols <- rainbow(length(my.taxas))
	ind <- grep("Actinobacteriota", my.phylums)
	if(length(ind)>0) my.cols[ind] <- colorRampPalette(c("mediumpurple4", "mediumorchid1"))(length(ind))
	ind <- grep("Bacteroidota", my.phylums)
	if(length(ind)>0) my.cols[ind] <- colorRampPalette(c("darkred", "chocolate1"))(length(ind))
	ind <- grep("Firmicutes", my.phylums)
	if(length(ind)>0) my.cols[ind] <- colorRampPalette(c("blue4", "cadetblue1"))(length(ind))
	ind <- grep("Proteobacteria", my.phylums)
	if(length(ind)>0) my.cols[ind] <- colorRampPalette(c("darkgreen", "darkseagreen1"))(length(ind))
	ind <- grep("Verrucomicrobiota", my.phylums)
	if(length(ind)>0) my.cols[ind] <- colorRampPalette(c("goldenrod", "khaki1"))(length(ind))
	my.unique.phy <- setdiff(my.unique.phy, c("Actinobacteriota", "Bacteroidota", "Firmicutes", "Proteobacteria", "Verrucomicrobiota")) 
	if(length(my.unique.phy)>0){
		my.col.lens <- c()
		my.col.ind <- c()
		for(i in 1:length(my.unique.phy)){
			temp.ind <- grep(my.unique.phy[i], my.phylums)
			my.col.ind <- c(my.col.ind, temp.ind)
			my.col.lens <- c(my.col.lens, length(temp.ind))
		}
		j <- 1
		k <- 1
		if(2*length(my.unique.phy)<=length(my.extreme.cols1)){
			for(i in 1:length(my.col.lens)){
        			my.cols[my.col.ind[k]:(my.col.ind[k]+my.col.lens[i]-1)] <- colorRampPalette(my.extreme.cols1[j:(j+1)])(my.col.lens[i])
				j <- j+2
				k <- k+my.col.lens[i]
			}
		} else stop("add colors in my.extreme.cols1")
	}

	my.phylums <- paste(my.phylums, my.taxas, sep="|")
	names(my.cols) <- my.phylums
	
	if(splitLegend==TRUE) layout(matrix(c(1,1,1,1,1,2,2,2,3,3,3), 1, 11)) else layout(matrix(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2), 1, 14))
	par(mar=c(mybottom, 1, 1, 0))

	if(is.element("other", colnames(ta.dfm.pct))){
		barplot(t(cbind(other.pct, my.ta.dfm.pct)), col=c("grey", my.cols), cex.names=0.7, cex.axis=0.8, border=NA, las=1,
			mgp=c(3,0.05,0), tcl=-0.05, main=mymain)
		mylegend <- c("other", my.phylums); myfill <- c("grey", my.cols[my.phylums])
	} else {
		barplot(t(my.ta.dfm.pct), col=my.cols, cex.names=0.7, cex.axis=0.8, border=NA, las=1,
			mgp=c(3,0.05,0), tcl=-0.05, main=mymain)
		mylegend <- my.phylums; myfill <- my.cols[my.phylums]
	}
	mylegend <- mylegend[length(mylegend):1]; myfill <- myfill[length(myfill):1]
	m <- round(length(mylegend)/2)
	if(splitLegend==TRUE){
		plot(1,1, type="n", ylim=c(0,100), xlim=c(0,50), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
		legend(x=0, y=100, legend=mylegend[1:m], fill=myfill[1:m], bty="n", border="white", cex=0.6)
		plot(1,1, type="n", ylim=c(0,100), xlim=c(0,50), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
		legend(x=0, y=100, legend=mylegend[(m+1):length(mylegend)], fill=myfill[(m+1):length(mylegend)], bty="n", border="white", cex=0.6)
	} else {
		plot(1,1, type="n", ylim=c(0,100), xlim=c(0,50), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
		legend(x=0, y=100, legend=mylegend, fill=myfill, bty="n", border="white", cex=0.6)
	}
}

get.splits <- function(ALL.classes, splits=100, prop.test=0.5){
	yes.group <- levels(ALL.classes)[1]

	y <- rep(0, length(ALL.classes))
	y[ALL.classes==yes.group] <- 1
	y <- as.integer(y)
	ind.0 <- which(y==0)
	ind.1 <- which(y==1)
	ind.test.m <- c()
	for(i in 1:splits){
		ind.test <- c(sample(ind.0, round(length(ind.0)*prop.test)), sample(ind.1, round(length(ind.1)*prop.test)))
		ind.test.m <- rbind(ind.test.m, ind.test)
	}
	return(ind.test.m)
}

#####################################################################
###### RANDOM FOREST
get.splits <- function(ALL.classes, splits=100, prop.test=0.5){
	yes.group <- levels(ALL.classes)[1]

	y <- rep(0, length(ALL.classes))
	y[ALL.classes==yes.group] <- 1
	y <- as.integer(y)
	ind.0 <- which(y==0)
	ind.1 <- which(y==1)
	ind.test.m <- c()
	for(i in 1:splits){
		ind.test <- c(sample(ind.0, round(length(ind.0)*prop.test)), sample(ind.1, round(length(ind.1)*prop.test)))
		ind.test.m <- rbind(ind.test.m, ind.test)
	}
	return(ind.test.m)
}

assess.categ <- function(pred.test, prob.cut.positive=0.5, positive.level="YES", ALL, ALL.classes){
	true.test <- ALL.classes; names(true.test) <- colnames(ALL); true.test <- true.test[names(pred.test)]
	pred.test.positive <- names(pred.test)[which(pred.test>prob.cut.positive)]
	pred.test.negative <- names(pred.test)[which(pred.test<=prob.cut.positive)]
	true.test.positive <- names(true.test)[which(true.test==positive.level)]
	true.test.negative <- names(true.test)[which(true.test!=positive.level)]
	VP <- length(intersect(pred.test.positive, true.test.positive))
	FP <- length(intersect(pred.test.positive, true.test.negative))
	VN <- length(intersect(pred.test.negative, true.test.negative))
	FN <- length(intersect(pred.test.negative, true.test.positive))
	res <- c(VP/(VP+FN), FP/(FP+VN)); names(res) <- c("VPR", "FPR")
	return(res)
}

assess.RF <- function(pred.test, prob.cut.positive=0.5, positive.level="YES", ALL, ALL.classes){
	true.test <- ALL.classes; names(true.test) <- colnames(ALL); true.test <- true.test[names(pred.test)]
	pred.test.positive <- names(pred.test)[which(pred.test>prob.cut.positive)]
	pred.test.negative <- names(pred.test)[which(pred.test<=prob.cut.positive)]
	true.test.positive <- names(true.test)[which(true.test==positive.level)]
	true.test.negative <- names(true.test)[which(true.test!=positive.level)]
	VP <- length(intersect(pred.test.positive, true.test.positive))
	FP <- length(intersect(pred.test.positive, true.test.negative))
	VN <- length(intersect(pred.test.negative, true.test.negative))
	FN <- length(intersect(pred.test.negative, true.test.positive))
	res <- c(VP/(VP+FN), FP/(FP+VN)); names(res) <- c("VPR", "FPR")
	return(res)
}

assessRFCateg <- function(ALL, ALL.classes, ind.test.m, ntree=1000, meta.lab="RF"){

	yes.group <- levels(ALL.classes)[1]
	mydata <- as.data.frame(t(ALL))
	mydata$y <- ALL.classes
	splits <- nrow(ind.test.m)
	my.mtry <- round(sqrt(ncol(mydata)-1))

	# calculamos los RF para cada split para las variables de mysel
	rf.pr.l <- vector("list", length=splits)
	RF.l <- vector("list", length=splits)
	for(i in 1:splits){
		inTrain <- setdiff(1:nrow(mydata), ind.test.m[i, ])
		training <- mydata[inTrain, ]
		testing <- mydata[-inTrain, ]
		RF <- randomForest(y ~ ., data=training, mtry=my.mtry, ntree=ntree, keep.forest=TRUE, importance=FALSE,
			xtest=testing[, setdiff(colnames(testing), "y")], ytest=testing[, "y"])
		# (si aquí la función randomForest da un error extraño, tenemos que asegurarnos de que los nombres de variables sean suficientemente simples
		# sin carácteres extraños)
		rf.pr.l[[i]] <- predict(RF, type="prob", newdata=testing)[, yes.group]
		RF.l[[i]] <- RF
		#print(i)
	}

	# curva ROC
	# Cortes de discriminación para el cálculo de la curva ROC
	cuts <- seq(from=0, to=1, length.out=50)
	VPR.v <- c(); FPR.v <- c()
	# evaluamos los RF para cada split y para cada cut
	for(mycut in cuts){
		vpr <- c(); fpr <- c()
		for(i in 1:splits){
			res <- assess.categ(pred.test=rf.pr.l[[i]], prob.cut.positive=mycut, positive.level=yes.group, ALL=ALL, ALL.classes=ALL.classes)
			vpr <- c(vpr, res["VPR"]); fpr <- c(fpr, res["FPR"])
		}
		VPR.v <- c(VPR.v, mean(vpr)); FPR.v <- c(FPR.v, mean(fpr))
		#print(mycut)
	}
	cuts <- format(cuts, digits=1)
	AUC <- abs(trapz(FPR.v, VPR.v))
	if(is.nan(AUC)){
		break()
	} else {
		AUC <- format(AUC, digits=2)
	}

	write.table(data.frame(FPR=FPR.v, VPR=VPR.v, discrimination.threshold=cuts), file=paste(meta.lab, "points.ROC", ncol(mydata)-1, "variables", "tsv", sep="."),
		row.names=FALSE, sep="\t", quote=FALSE)
	mymain <- paste(splits, "splits ", ncol(mydata)-1, "variables ", "AUC:", AUC, sep=" ")
	pdf(file=paste(meta.lab, "ROC", ncol(mydata)-1, "variables", "pdf", sep="."))
	plot(x=FPR.v, y=VPR.v, type="n", xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), main=mymain)
	abline(a=0, b=1, col="red")
	text(x=FPR.v, y=VPR.v, labels=cuts, col="blue", cex=0.5)
	lines(x=FPR.v, y=VPR.v, col="blue")
	dev.off()

	return(AUC)
	# A modo de guía para interpretar las curvas ROC se han establecido los siguientes intervalos para los valores de AUC:
	# [0.5, 0.6): Test malo.
	# [0.6, 0.75): Test regular.
	# [0.75, 0.9): Test bueno.
	# [0.9, 0.97): Test muy bueno.
	# [0.97, 1): Test excelente.
}

doRF <- function(ALL, ALL.classes, meta.label="", description=NULL, splits=100, prop.test=0.5, ntree=1000, quantile.start=0.9, num.confirmed=40){

#################
## Lectura de datos:
## ALL: matriz con features en rows y samples en columns (ha de estar normalizada por samples: por ejemplo,
##		la suma total aplicada a cada sample ha de ser aprox. constante)
## ALL.classes: factor (variable categórica) con la clasificación de las samples (mismo orden que las columnas de ALL)
#load(file="ALL.RData")
#load(file="ALL.classes.RData")

rownames(ALL) <- sub("/", ".", rownames(ALL))

yes.group <- levels(ALL.classes)[1]

# Conviene tener filt==TRUE si el número de predictores candidatos originales es > 800 (se eliminarán variables con baja variabilidad)
#filt <- TRUE
if(nrow(ALL)>800) filt <- TRUE else filt <- FALSE

#################
## (1)
## Construcción del cross-validation design
# fijamos el número de particiones training/test (número de splits):
#splits <- 100
# proporción en cada grupo de muestras test
#prop.test <- 0.5

#if(length(system("ls ind.test.m.RData", intern=TRUE))>0){
#	load(file="ind.test.m.RData")
#} else {
	y <- rep(0, length(ALL.classes))
	y[ALL.classes==yes.group] <- 1
	y <- as.integer(y)
	ind.0 <- which(y==0)
	ind.1 <- which(y==1)
	ind.test.m <- c()
	for(i in 1:splits){
		ind.test <- c(sample(ind.0, round(length(ind.0)*prop.test)), sample(ind.1, round(length(ind.1)*prop.test)))
		ind.test.m <- rbind(ind.test.m, ind.test)
	}
	#save(ind.test.m, file="ind.test.m.RData")
#}

#################
## (2) y (3)
## RF

# Conviene tener filt==TRUE si el número de predictores candidatos originales es > 800
# Nos quedamos con el 60% con mayor variabilidad (asumimos que las variables de entrada ya han pasado
# un filtro por señal que garantice que cada una de ellas tiene algo de señal decente en al menos uno de los grupos de samples).
#if(filt==TRUE){
#	cond <- apply(ALL, MARGIN=1, sd)
#	myvars <- rownames(ALL)[which(cond>quantile(cond, probs=0.2))]
#} else {
	myvars <- rownames(ALL)
#}

mydata <- as.data.frame(t(ALL[myvars, ]))
mydata$y <- ALL.classes

splits <- nrow(ind.test.m)

my.mtry <- round(sqrt(ncol(mydata)-1))

# calculamos los RF para cada split 
rf.pr.l <- vector("list", length=splits)
RF.l <- vector("list", length=splits)
for(i in 1:splits){
	inTrain <- setdiff(1:nrow(mydata), ind.test.m[i, ])
	training <- mydata[inTrain, ]
	testing <- mydata[-inTrain, ]
	RF <- randomForest(y ~ ., data=training, mtry=my.mtry, ntree=ntree, keep.forest=TRUE, importance=TRUE,
		xtest=testing[, setdiff(colnames(testing), "y")], ytest=testing[, "y"])
	RF.l[[i]] <- RF
	rf.pr.l[[i]] <- predict(RF, type="prob", newdata=testing)[, yes.group]
	print(i)
}
# (si aquí la función randomForest da un error extraño, tenemos que asegurarnos de que los nombres de variables sean suficientemente simples
# sin carácteres extraños)

# Cortes de discriminación para el cálculo de la curva ROC
cuts <- seq(from=0, to=1, length.out=50)
VPR.v <- c(); FPR.v <- c()
# evaluamos los RF para cada split y para cada cut
for(mycut in cuts){
	vpr <- c(); fpr <- c()
	for(i in 1:splits){
		res <- assess.RF(pred.test=rf.pr.l[[i]], prob.cut.positive=mycut, positive.level=yes.group, ALL=ALL, ALL.classes=ALL.classes)
		vpr <- c(vpr, res["VPR"]); fpr <- c(fpr, res["FPR"])
	}
	VPR.v <- c(VPR.v, mean(vpr)); FPR.v <- c(FPR.v, mean(fpr))
	print(mycut)
}

cuts <- format(cuts, digits=1)
AUC <- abs(trapz(FPR.v, VPR.v))

if(is.nan(AUC)){
	warning("El AUC sale NaN")
	return()
} else {
	AUC <- format(AUC, digits=2)
}
#write.table(data.frame(FPR=FPR.v, VPR=VPR.v), file=paste("points.ROC.RF", ncol(mydata)-1, "variables", meta.label, "tsv", sep="."), row.names=FALSE, sep="\t", quote=FALSE)
mymain <- paste(splits, "splits ", ncol(mydata)-1, "variables ", "AUC:", AUC, sep=" ")
#pdf(file=paste("ROC.RF", ncol(mydata)-1, "variables", meta.label, "pdf", sep="."))
#plot(x=FPR.v, y=VPR.v, type="n", xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), main=mymain)
#abline(a=0, b=1, col="red")
#text(x=FPR.v, y=VPR.v, labels=cuts, col="blue", cex=0.6)
#lines(x=FPR.v, y=VPR.v, col="blue")
#dev.off()

# Cálculo de importance
RF <- RF.l[[1]]
impor <- importance(RF)
impor.acc.m <- impor[, "MeanDecreaseAccuracy"]
impor.gini.m <- impor[, "MeanDecreaseGini"]
# Most significant variables in the model
for(i in 2:splits){
RF <- RF.l[[i]]
impor <- importance(RF)
impor.acc.v <- impor[, "MeanDecreaseAccuracy"]
impor.gini.v <- impor[, "MeanDecreaseGini"]
impor.acc.m <- cbind(impor.acc.m, impor.acc.v)
impor.gini.m <- cbind(impor.gini.m, impor.gini.v)
}
rownames(impor.acc.m) <- rownames(ALL)
rownames(impor.gini.m) <- rownames(ALL)
impor.acc <- apply(impor.acc.m, MARGIN=1, mean)
impor.gini <- apply(impor.gini.m, MARGIN=1, mean)
impor <- data.frame(MeanDecreaseAccuracy=impor.acc, MeanDecreaseGini=impor.gini)

mytb <- impor[order(impor[, "MeanDecreaseAccuracy"], decreasing=TRUE), c("MeanDecreaseAccuracy", "MeanDecreaseGini")]
if(!is.null(description)) mytb$description <- description[rownames(mytb)]
write.csv2(mytb, file=paste("importance.MDA", meta.label, "csv", sep="."))
mytb$feature <- rownames(mytb)
mytb$decision <- rep("Rejected", nrow(mytb))
mytb$decision[1:num.confirmed] <- "Confirmed"
mytb <- mytb[, c("feature", "decision", "MeanDecreaseAccuracy", "MeanDecreaseGini", "description")]
write.table(mytb, file=paste("importance.MDA", meta.label, "tsv", sep="."), sep="\t", quote=FALSE, row.names=FALSE)
mytb <- impor[order(impor[, "MeanDecreaseGini"], decreasing=TRUE), c("MeanDecreaseGini", "MeanDecreaseAccuracy")]
if(!is.null(description)) mytb$description <- description[rownames(mytb)]
#write.csv2(mytb, file=paste("importance.MDG", meta.label, "csv", sep="."))

impor <- impor[order(impor[, "MeanDecreaseAccuracy"], decreasing=TRUE), ]

## Obtener distintos sets de variables con número decreciente de variables (definiendo los cortes de probs)
# probs: 0.95, 0.97, 0.98, 0.99, 0.995
cut1 <- quantile(impor[, "MeanDecreaseAccuracy"], probs=quantile.start)
#cut2 <- quantile(impor[, "MeanDecreaseGini"], probs=quantile.start)
#mysel <- intersect(rownames(impor)[which(impor[, "MeanDecreaseAccuracy"]>=cut1)], rownames(impor)[which(impor[, "MeanDecreaseGini"]>=cut2)])
mysel <- rownames(impor)[which(impor[, "MeanDecreaseAccuracy"]>=cut1)]
#write(mysel, file=paste("RF", length(mysel), "important.variables", meta.label, "txt", sep="."), sep="\n")

# A modo de guía para interpretar las curvas ROC se han establecido los siguientes intervalos para los valores de AUC:
# [0.5, 0.6): Test malo.
# [0.6, 0.75): Test regular.
# [0.75, 0.9): Test bueno.
# [0.9, 0.97): Test muy bueno.
# [0.97, 1): Test excelente.

## Repetir esta parte para distintos mysel hasta que no haya ganancia de curva ROC
# Evaluamos el modelo a partir de la selección hecha a partir del primer cálculo de importance

#myseq <- seq(from=quantile.start, to=0.995, length.out=quantile.length.out)
#myseq <- myseq[2:length(myseq)]
myseq <- seq(from=length(mysel)-1, to=2, by=-1)
AUC.old <- as.numeric(AUC)
AUC.new <- AUC.old
k <- 1

res.df <- c()

#while(AUC.new>=AUC.old&length(mysel)>2&k<=length(myseq)){
while(k<=length(myseq)){
myvars <- mysel
mydata <- as.data.frame(t(ALL[myvars, ]))
mydata$y <- ALL.classes
my.mtry <- round(sqrt(ncol(mydata)-1))
# calculamos los RF para cada split 
rf.pr.l <- vector("list", length=splits)
RF.l <- vector("list", length=splits)
for(i in 1:splits){
	inTrain <- setdiff(1:nrow(mydata), ind.test.m[i, ])
	training <- mydata[inTrain, ]
	testing <- mydata[-inTrain, ]
	RF <-randomForest(y ~ ., data=training, mtry=my.mtry, ntree=ntree, keep.forest=TRUE, importance=FALSE,
		xtest=testing[, setdiff(colnames(testing), "y")], ytest=testing[, "y"])
	rf.pr.l[[i]] <- predict(RF, type="prob", newdata=testing)[, yes.group]
	RF.l[[i]] <- RF
	print(i)
}
cuts <- seq(from=0, to=1, length.out=50)
VPR.v <- c(); FPR.v <- c()
# evaluamos los RF para cada split y para cada cut
for(mycut in cuts){
	vpr <- c(); fpr <- c()
	for(i in 1:splits){
		res <- assess.RF(pred.test=rf.pr.l[[i]], prob.cut.positive=mycut, positive.level=yes.group, ALL=ALL, ALL.classes=ALL.classes)
		vpr <- c(vpr, res["VPR"]); fpr <- c(fpr, res["FPR"])
	}
	VPR.v <- c(VPR.v, mean(vpr)); FPR.v <- c(FPR.v, mean(fpr))
	print(mycut)
}
#cuts <- format(cuts, digits=1)
AUC <- abs(trapz(FPR.v, VPR.v))
if(is.nan(AUC)){
	break()
} else {
	AUC <- format(AUC, digits=2)
}
if(length(mysel)<30){
#write.table(data.frame(FPR=FPR.v, VPR=VPR.v), file=paste("points.ROC.RF", ncol(mydata)-1, "variables", meta.label, "tsv", sep="."), row.names=FALSE, sep="\t", quote=FALSE)
mymain <- paste(splits, "splits ", ncol(mydata)-1, "variables ", "AUC:", AUC, sep=" ")
#pdf(file=paste("ROC.RF", ncol(mydata)-1, "variables", meta.label, "pdf", sep="."))
#plot(x=FPR.v, y=VPR.v, type="n", xlab="FPR", ylab="TPR", xlim=c(0,1), ylim=c(0,1), main=mymain)
#abline(a=0, b=1, col="red")
#text(x=FPR.v, y=VPR.v, labels=cuts, col="blue", cex=0.6)
#lines(x=FPR.v, y=VPR.v, col="blue")
#dev.off()
}
AUC.old <- AUC.new
#cut1 <- quantile(impor[, "MeanDecreaseAccuracy"], probs=myseq[k])
#cut2 <- quantile(impor[, "MeanDecreaseGini"], probs=myseq[k])
#mysel <- intersect(rownames(impor)[which(impor[, "MeanDecreaseAccuracy"]>=cut1)], rownames(impor)[which(impor[, "MeanDecreaseGini"]>=cut2)])
mysel <- rownames(impor)[1:myseq[k]]
#if(length(mysel)<30) write(mysel, file=paste("RF", length(mysel), "important.variables", meta.label, "txt", sep="."), sep="\n")
print(k)
AUC.new <- as.numeric(AUC)
k <- k+1
res.df <- rbind(res.df, data.frame(AUC=as.numeric(AUC), n.vars=ncol(mydata)-1, file.vars=paste("RF", ncol(mydata)-1, "important.variables", meta.label, "txt", sep=".")))
}

write.table(res.df, file=paste("results.RF", meta.label, "tsv", sep="."), sep="\t", row.names=FALSE, quote=FALSE)

pdf(file=paste("results.RF", meta.label, "pdf", sep="."))
plot(x=res.df$n.vars, y=res.df$AUC, pch=19, col="blue", xlim=c(max(res.df$n.vars), min(res.df$n.vars)))
lines(x=res.df$n.vars, y=res.df$AUC, col="blue")
abline(v=10, col="red")
dev.off()

return()
################################
# Qué haríamos para predecir el estatus de nuevas muestras basándose en un test que observara sólo unos pocos biomarcadores, digamos ¿15?

# (1) Se vuelven a entrenar los 50 RF usando sólo los biomarcadores seleccionados y a partir de unos nuevos 50 splits con proporción
#	mayor para training (se puede hacer un leave-one-out) (ya no necesitamos test sets, ya que asumimos como accuracy del modelo
#	la obtenida anteriormente para esos 15 predictores).
# (2) Para cada RF se usa como test la nueva (o nuevas) variables a clasificar
# (3) Como clasificación final se coge la media de las 50 clasificaciones

}

#####################################################################
###### OPLS-regression cross-validation
assess.OPLS.regression <- function(M, Y, IV, len.train, splits){
	res <- c()
	for(i in 1:splits){
		train <- sample(1:ncol(M), len.train)
		opls.res <- opls(x=t(M[IV, ]), y=Y, predI=1, orthoI=1, subset=train, plotL=FALSE)
		# devolvemos el error relativo ||v-v*|| / ||v|| donde v* es la estimación del verdadero v
		new.res <- dist(rbind(Y[-train], predict(opls.res, t(M[IV, ])[-train, ])))
		new.res <- new.res/dist(rbind(Y[-train], rep(0, length(Y)-length(train))))
		res <- c(res, new.res)
	}
	return(median(res))
}

select.vars.OPLS.regression.quantile <- function(M, Y, vip, splits, len.train, quantile.start, quantile.by, error.ini=1000){

	mysel <- names(vip)[which(vip>=quantile(vip, probs=quantile.start))]
	
	myseq <- seq(from=quantile.start, to=0.995, by=quantile.by)
	myseq <- myseq[2:length(myseq)]
	error.old <- error.ini
	error.new <- error.old
	k <- 1

	res.df <- c()
	probs <- quantile.start

	while(error.new<=error.old&length(mysel)>2&k<=length(myseq)){
		error <- assess.OPLS.regression(M=M, Y=Y, IV=mysel, len.train=len.train, splits=splits)
		error.old <- error.new
		res.df <- rbind(res.df, data.frame(error=error, n.vars=length(mysel), probs=probs, vars=paste(mysel, collapse=","), stringsAsFactors=FALSE))
		mysel <- names(vip)[which(vip>=quantile(vip, probs=myseq[k]))]
		error.new <- error
		probs <- myseq[k]
		k <- k+1
	}
	
	return(res.df)
}

select.vars.OPLS.regression <- function(M, Y, vip, splits, len.train, by, error.ini=1000){

	mysel <- names(vip)
	
	myseq <- seq(from=length(vip), to=3, by=by)
	myseq <- myseq[2:length(myseq)]
	error.old <- error.ini
	error.new <- error.old
	k <- 1

	res.df <- c()

	while(error.new<=error.old&length(mysel)>2&k<=length(myseq)){
		error <- assess.OPLS.regression(M=M, Y=Y, IV=mysel, len.train=len.train, splits=splits)
		error.old <- error.new
		res.df <- rbind(res.df, data.frame(error=error, n.vars=length(mysel), vars=paste(mysel, collapse=","), stringsAsFactors=FALSE))
		mysel <- names(vip)[1:myseq[k]]
		error.new <- error
		k <- k+1
	}
	
	return(res.df)
}


#####################################################################
###### randomForest con respuesta continua cross-validation
assess.randomForest.regression <- function(M, Y, IV, len.train, splits, ntree, importance=FALSE){

	mydata <- as.data.frame(t(M[IV, ]))
	mydata$y <- Y
	my.mtry <- round(sqrt(ncol(mydata)-1))
	impor.acc.m <- c(); impor.gini.m <- c()
	res <- c()

	for(i in 1:splits){
		inTrain <- sample(1:nrow(mydata), len.train)
		training <- mydata[inTrain, ]
		testing <- mydata[-inTrain, ]
		RF <- randomForest(y ~ ., data=training, mtry=my.mtry, ntree=ntree, keep.forest=TRUE, importance=importance,
						xtest=testing[, setdiff(colnames(testing), "y")], ytest=testing[, "y"])
		new.res <- dist(rbind(Y[-inTrain], predict(RF, newdata=testing)))
		new.res <- new.res/dist(rbind(Y[-inTrain], rep(0, length(Y)-length(inTrain))))
		res <- c(res, new.res)
		if(importance==TRUE){
			impor <- importance(RF)
			impor.acc.v <- impor[, 1]
			impor.gini.v <- impor[, 2]
			impor.acc.m <- cbind(impor.acc.m, impor.acc.v)
			impor.gini.m <- cbind(impor.gini.m, impor.gini.v)
		}
		print(i)
	}
	
	if(importance==TRUE){
		rownames(impor.acc.m) <- IV
		rownames(impor.gini.m) <- IV
		impor.acc <- apply(impor.acc.m, MARGIN=1, FUN=function(x) median(x, na.rm=TRUE))
		impor.gini <- apply(impor.gini.m, MARGIN=1, FUN=function(x) median(x, na.rm=TRUE))
		impor <- data.frame(MeanDecreaseAccuracy=impor.acc, MeanDecreaseGini=impor.gini)
		mytb <- impor[order(impor[, "MeanDecreaseAccuracy"], decreasing=TRUE), c("MeanDecreaseAccuracy", "MeanDecreaseGini")]
		mytb$DESCRIPTION <- DESCRIPTION[rownames(mytb)]
		mytb.MDA <- mytb
		mytb <- impor[order(impor[, "MeanDecreaseGini"], decreasing=TRUE), c("MeanDecreaseGini", "MeanDecreaseAccuracy")]
		mytb$DESCRIPTION <- DESCRIPTION[rownames(mytb)]
		mytb.MDG <- mytb
	} else {
		mytb.MDA <- mytb.MDG <- NULL
	}

	return(list(error=median(res), mytb.MDA=mytb.MDA, mytb.MDG=mytb.MDG))
}

select.vars.randomForest.regression <- function(M, Y, vip.MDG, vip.MDA, splits, len.train, quantile.start, quantile.by, error.ini, ntree){

	mysel.MDG <- names(vip.MDG)[which(vip.MDG>=quantile(vip.MDG, probs=quantile.start))]
	mysel.MDA <- names(vip.MDA)[which(vip.MDA>=quantile(vip.MDA, probs=quantile.start))]
	mysel <- intersect(mysel.MDG, mysel.MDA)
	
	myseq <- seq(from=quantile.start, to=0.995, by=quantile.by)
	myseq <- myseq[2:length(myseq)]
	error.old <- error.ini
	error.new <- error.old
	k <- 1

	res.df <- c()
	probs <- quantile.start

	while(error.new<=error.old&length(mysel)>2&k<=length(myseq)){
		error <- assess.randomForest.regression(M=M, Y=Y, IV=mysel, len.train=len.train, splits=splits, ntree=ntree)$error
		error.old <- error.new
		res.df <- rbind(res.df, data.frame(error=error, n.vars=length(mysel), probs=probs, vars=paste(mysel, collapse=","), stringsAsFactors=FALSE))
		mysel.MDG <- names(vip.MDG)[which(vip.MDG>=quantile(vip.MDG, probs=myseq[k]))]
		mysel.MDA <- names(vip.MDA)[which(vip.MDA>=quantile(vip.MDA, probs=myseq[k]))]
		mysel <- intersect(mysel.MDG, mysel.MDA)
		error.new <- error
		probs <- myseq[k]
		k <- k+1
	}
	
	return(res.df)
}

# NOTAS sobre opls:
# x: Numerical data frame or matrix (observations x variables)
# y: a factor (same length as 'x' row number) for single response
# predI: Integer: number of components (predictive componenents in
#          case of PLS and OPLS) to extract; for OPLS, predI is
#          (automatically) set to 1
#  orthoI: Integer: number of orthogonal components (for OPLS only);
#          when set to 0 [default], PLS will be performed; otherwise
#          OPLS will be peformed; when set to NA, OPLS is performed and
#          the number of orthogonal components is automatically computed
#          by using the cross-validation (with a maximum of 9 orthogonal
#          components).
# Error bastante habitual en opls:
#Error: No model was built because the first orthogonal component was already not significant;
#Select a number of orthogonal components of 1 if you want the algorithm to compute a model despite this.

assessOPLSCateg <- function(ALL, ALL.classes, ind.test.m){

	yes.group <- levels(ALL.classes)[1]
	no.group <- levels(ALL.classes)[2]
	splits <- nrow(ind.test.m)

	# calculamos los OPLS para cada split para las variables de mysel
	OPLS.pr.l <- vector("list", length=splits)
	OPLS.l <- vector("list", length=splits)
	for(i in 1:splits){
		inTrain <- setdiff(1:ncol(ALL), ind.test.m[i, ])
		OPLS <- opls(x=t(ALL), y=ALL.classes, predI=1, orthoI=1, subset=inTrain, plotL=FALSE)
		OPLS.pr.l[[i]] <- predict(OPLS, t(ALL)[-inTrain, ])
		OPLS.l[[i]] <- OPLS
		#print(i)
	}

	# predict.opls no tiene implementado el type="prob", así que no podemos representar la curva ROC
	# (tampoco sabemos el corte de discriminación que usan en su método)
	# simplemente, representaremos la distribución de VPR y FDR dada por corte de discriminación usado por predict.opls
	VPR <- c(); FPR <- c()
	# evaluamos los OPLS para cada split y para cada cut
	for(i in 1:splits){
		tem <- OPLS.pr.l[[i]]
		YES <- rep(NA, length(tem))
		YES[which(tem==no.group)] <- 0.4; YES[which(tem==yes.group)] <- 0.6
		names(YES) <- names(tem)
		res <- assess.categ(pred.test=YES, prob.cut.positive=0.5, positive.level=yes.group, ALL=ALL, ALL.classes=ALL.classes)
		VPR <- c(VPR, res["VPR"]); FPR <- c(FPR, res["FPR"])
	}
	
	return(list(VPR=VPR, FPR=FPR))
}

assessOPLSregression <- function(M, Y, IV, ind.test.m){
	splits <- nrow(ind.test.m)
	res <- c()
	for(i in 1:splits){
		train <- setdiff(1:nrow(mydata), ind.test.m[i, ])
		opls.res <- opls(x=t(M[IV, ]), y=Y, predI=1, orthoI=1, subset=train, plotL=FALSE)
		# devolvemos el error relativo ||v-v*|| / ||v|| donde v* es la estimación del verdadero v
		new.res <- dist(rbind(Y[-train], predict(opls.res, t(M[IV, ])[-train, ])))
		new.res <- new.res/dist(rbind(Y[-train], rep(0, length(Y)-length(train))))
		res <- c(res, new.res)
	}
	return(res)
}

assessRFregression <- function(M, Y, IV, ind.test.m, ntree){

	splits <- nrow(ind.test.m)
	mydata <- as.data.frame(t(M[IV, ]))
	mydata$y <- Y
	my.mtry <- round(sqrt(ncol(mydata)-1))
	impor.acc.m <- c(); impor.gini.m <- c()
	res <- c()

	for(i in 1:splits){
		inTrain <- setdiff(1:nrow(mydata), ind.test.m[i, ])
		training <- mydata[inTrain, ]
		testing <- mydata[-inTrain, ]
		RF <- randomForest(y ~ ., data=training, mtry=my.mtry, ntree=ntree, keep.forest=TRUE, importance=FALSE,
						xtest=testing[, setdiff(colnames(testing), "y")], ytest=testing[, "y"])
		# devolvemos el error relativo ||v-v*|| / ||v|| donde v* es la estimación del verdadero v
		new.res <- dist(rbind(Y[-inTrain], predict(RF, newdata=testing)))
		new.res <- new.res/dist(rbind(Y[-inTrain], rep(0, length(Y)-length(inTrain))))
		res <- c(res, new.res)
		#print(i)
	}

	return(res)
}


###########################################################################################################
###########################################################################################################
###### ANCOM-II
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5682008/
## "Analysis of Microbiome Data in the Presence of Excess Zeros"
## Metodología para detectar features con abundancia diferencial en 2 o más grupos en datos composicionales.
## Está pensado en origen para datos composicionales: para tener en cuenta la estructura simplex inherente a este tipo de datos, realizan una transformación
## basada en log-ratios. Como ventaja frente a otras aproximaciones, no hace asunciones sobre la distribución que siguen los datos. En este sentido puede
## resultar más apropiado que ALDEX2, aunque, se le parece bastante en lo relativo a la transformación CLR de los datos y al tratamiento de los ceros outliers.
## Esta metodología supera, tanto en FDR como en falsos negativos a la tradicional t-test o wilcox-test sobre los datos normalizados en porcentajes.
## Por otro lado, la técnica de sumar 1's (pseudocounts) seguido de la tradicional t-test o wilcox-test sobre los datos normalizados en porcentajes, supera
## a ANCOM-II en FDR, pero arroja muchos más falsos negativos.
## En cuanto a la comparación de ANCOM-II frente a deseq2, son iguales en FDR pero ANCOM-II supera deseq2 en falsos negativos. Además deseq2, aunque usa
## una distribución apropiada para modelar datos composicionales y tiene la ventaja de incorporar métodos de inferencia bayesiana, no tiene en cuenta la
## estructura simplex inherente a este tipo de datos.
## Por todo ello, se recomienda escoger ANCOM-II frente a deseq2 para la selección de features taxonómicas diferencialmente abundantes.
## Tener en cuenta, por otro lado que deseq2 está testado en datos de rna-seq, donde el número de features supera bastante al caso de las matrices del 16S, por tanto,
## como los sesgos derivados de la imposición de la estructura simplex en composiciones disminuyen a medida que aumenta el número de features, es posible que en contextos
## con elevado número de features deseq2 sea equivalente a ANCOM-II.

#A.	Clustering algorithm
pop_detect <- function(obs){

  TOL = 1e-2;
  MAXITS=1000
  its = 1;
  change = 1;
  obs=obs
  q1=as.numeric(quantile(obs, 0.25))
  q3=as.numeric(quantile(obs, 0.75))
  iqr=q3-q1
  mu1old=q1
  mu2old=q3
  si1old=iqr/2
  si2old=iqr/2
  yold=ifelse((obs<q1+(iqr/2)), 0,1)
  piold=length(which(yold==0))/length(obs)
  while((its <= MAXITS) && (change > TOL) ){
    ynew=ifelse((dnorm(obs, mean=mu1old,sd=si1old)/dnorm(obs, mean=mu2old,sd=si2old))>((1-piold)/piold), 0,1)
    obs_frame=data.frame(obs,ynew)
    obsnew1=filter(obs_frame, obs_frame$ynew==0)[,1]
    obsnew2=filter(obs_frame, obs_frame$ynew==1)[,1]
    pinew=length(which(ynew==0))/length(obs)
    if(pinew>0.999){pinew=0.999}
    if(pinew<0.001){pinew=0.001}
    mu1new=ifelse(length(obsnew1)==0, mean(obsnew2), mean(obsnew1))
    si1new=ifelse(length(obsnew1)==1, 10^(-4),sd(obsnew1))
    si1new=ifelse(length(obsnew1)==0, sd(obsnew2), si1new)
    mu2new=ifelse(length(obsnew2)==0,mean(obsnew1),mean(obsnew2))
    si2new=ifelse(length(obsnew2)==1, 10^(-4), sd(obsnew2))
    si2new=ifelse(length(obsnew2)==0, sd(obsnew1), si2new)
    change = max(abs(mu1new-mu1old),abs(mu2new-mu2old),abs(pinew-piold),abs(si1old-si1new),abs(si2old-si2new)) ;
    mu1old=mu1new
    mu2old=mu2new
    si1old=si1new
    si2old=si2new
    piold=pinew
    its = its + 1;}
  if(its == MAXITS + 1){
    sprintf('max iterations');}
  return(list(mu1new,mu2new,si1new,si2new,pinew, obs, ynew, its));

}

####STEP A1
struc_zero <- function(OTUdat = otu, Vardat = meta, main.var = "population", p=NULL){

  if(is.null(p)==T){p=0.05}
main_var=main.var ;  meta=Vardat
ind=c(which(names(meta)=="Sample.ID"),which(names(meta)==main_var))
red_map=meta[,ind];   merged=merge(red_map, OTUdat, by="Sample.ID"); names(merged)[2]="pop"
#  merged_implement=merged[,-which(names(merged)=="Sample.ID")]
full_dat=merged;  datasplit<-list(); sid_list=list()
  for (i in 1:length(unique(full_dat$pop))){datasplit[[i]]=filter(full_dat, pop==unique(full_dat$pop)[i])
    sid_list[[i]]= datasplit[[i]]$Sample.ID; datasplit[[i]]=datasplit[[i]] %>% select(-Sample.ID)}
sid=unlist(sid_list); d=dim(datasplit[[1]])[2]; zeros=matrix(0,length(unique(full_dat$pop)),(d-1))
  for (i in 1:length(unique(full_dat$pop))){for (j in 2:d){
      microbe1=datasplit[[i]][,j];  microbe=na.omit(microbe1)
      x=sum(microbe!=0);    n=length(microbe) ;  r=x/n;      
      if(r<=p){datasplit[[i]][,j]=NA};  if(r<=p){zeros[i,(j-1)]=1}}}
dat_adjust_struc=datasplit[[1]]
  for (i in 2:length(unique(full_dat$pop))){dat_adjust_struc=rbind(dat_adjust_struc,datasplit[[i]])}
dat_adjust=data.frame(Sample.ID=sid,dat_adjust_struc)
rownames(zeros)=unique(full_dat$pop)
return(list(dat_adjust,zeros))

}

reset_dat_missing <- function(dat,pii=NULL, ref_name=NULL, zeros){

  if(is.null(ref_name)==T){normal_frame=normalizer_gm(dat)}
  if(is.null(ref_name)==F){normal_frame=normalizer_ref(dat,ref_name)}
  if(is.null(pii)==T){pii=0.25}; 
    #if(is.null(ref)==T){ref=nz_ref(dat)}
  ####ref_name is REFERENCE MICORBE###
  dat=merge(normal_frame,dat,by="Sample.ID")
  datasplit<-list();sid_list=list();pop_list=list(); normal_list=list()
  for (i in 1:length(unique(dat$pop))){datasplit[[i]]=filter(dat, pop==unique(dat$pop)[i])
  sid_list[[i]]= datasplit[[i]]$Sample.ID; pop_list[[i]]=datasplit[[i]]$pop
  normal_list[[i]]= datasplit[[i]]$normal; datasplit[[i]]=datasplit[[i]] %>% select(-c(pop,Sample.ID,normal))}
  popu=unlist(pop_list); sid=unlist(sid_list); normalu=unlist(normal_list)
  sub_dat=datasplit;   sub_dat_adjust=sub_dat
  d=dim(sub_dat[[1]])[2] ;  adj_microbes=matrix(0,length(unique(dat$pop)),d)
  
  for (i in 1:length(unique(dat$pop))){sub_dat_pseudo=sub_dat[[i]]+1 ;    nsub=dim(sub_dat_pseudo)[1]
  d=dim(sub_dat_pseudo)[2];     normalizer=matrix(0,nsub,1)
  normalizer=as.vector(normal_list[[i]])
  all_ind=which(colSums(zeros)==0);  
  if(is.null(ref_name)==T){ind_adjust=all_ind};  
  if(is.null(ref_name)==F){ind_adjust=setdiff(all_ind, which(names(sub_dat_pseudo)==ref_name))}
  for (j in ind_adjust){
  microbe=log(as.numeric(sub_dat_pseudo[,j]))-normalizer
  aa=pop_detect(microbe) ;    cluster_seq=aa[[7]]
  par=c(aa[[1]],aa[[2]],aa[[3]],aa[[4]],aa[[5]]) ;  rightend=par[[1]]+(1.96)*par[[3]]
  leftend=par[[2]]-(1.96)*par[[4]];     pileft=par[[5]];  piright=1-par[[5]];   rml=0;
  if(pileft<pii){rml=1} ;  
  if((rightend<leftend) && (rml==1)){if(length(cluster_seq==0)>0){ 
  adj_microbes[i,j]=1; clus_ind=which(cluster_seq==0) ; sub_dat_adjust[[i]][clus_ind,j]=NA}}}}
   
  number_adj_microbes=list() ;  names_adj_microbes=list()
  
  for (i in 1:length(unique(dat$pop))){number_adj_microbes[[i]]=sum(adj_microbes[i,])
  
  names_adj_microbes[[i]]=names(datasplit[[i]])[which(adj_microbes[i,]==1)]}
  
  dat_adjust=sub_dat_adjust[[1]]
  for (i in 2:length(unique(dat$pop))){dat_adjust=rbind(dat_adjust,sub_dat_adjust[[i]])}
  
  dat_adjust=data.frame(Sample.ID=sid,pop=popu, normal=normalu,dat_adjust)
  
  out=list(dat_adjust, number_adj_microbes, names_adj_microbes,adj_microbes)
  names(out)=c("adjusted data", "number of adjusted microbes", "names of adjusted microbes", "adjusted microbes")
  return(out)

}

normalizer_ref <- function(dat, ref_name){

  microbiome=dat %>% select(-c(Sample.ID,pop))
  nbig=dim(microbiome)[1];    d=dim(microbiome)[2]
  datasplit<-list();    sid_list=list()
  for (i in 1:length(unique(dat$pop))){ datasplit[[i]]=filter(dat, pop==unique(dat$pop)[i])
  sid_list[[i]]=datasplit[[i]]$Sample.ID;      datasplit[[i]]=datasplit[[i]] %>% select(-c(Sample.ID,pop))
  datasplit[[i]]=datasplit[[i]]+1};
  mean_ref_logscale=matrix(0,length(unique(dat$pop)),1);    normalize_logscale=list() 
  for (i in 1:length(unique(dat$pop))){
    mean_ref_logscale[i]=mean(log(as.matrix(datasplit[[i]][,ref_name])),na.rm=T)
    normalize_logscale[[i]]=log(as.matrix(datasplit[[i]][,ref_name]))-mean_ref_logscale[i]}
  normalize_ref=cbind(sid_list[[1]],normalize_logscale[[1]])
  for (i in 2:length(unique(dat$pop))){
    temp=cbind(sid_list[[i]],normalize_logscale[[i]]);      normalize_ref=rbind(normalize_ref,temp)}
  normalize_ref=data.frame(normalize_ref);   names(normalize_ref)=c("Sample.ID","normal")
  return(normalize_ref)

}

normalizer_gm <- function(dat){

    microbiome=dat %>% select(-c(Sample.ID,pop))
    nbig=dim(microbiome)[1];    d=dim(microbiome)[2]
    datasplit<-list();    sid_list=list()
    for (i in 1:length(unique(dat$pop))){ datasplit[[i]]=filter(dat, pop==unique(dat$pop)[i])
      sid_list[[i]]=datasplit[[i]]$Sample.ID;      datasplit[[i]]=datasplit[[i]] %>% select(-c(Sample.ID,pop))
      datasplit[[i]]=datasplit[[i]]+1};
    mean_ref_logscale=matrix(0,length(unique(dat$pop)),1);    normalize_logscale=list() 
  for (i in 1:length(unique(dat$pop))){
      mean_ref_logscale[i]=mean(rowMeans(log(as.matrix(datasplit[[i]])),na.rm=T),na.rm=T)
      normalize_logscale[[i]]=rowMeans(log(as.matrix(datasplit[[i]])),na.rm=T)-mean_ref_logscale[i]}
      normalize_gm=cbind(sid_list[[1]],normalize_logscale[[1]])
  for (i in 2:length(unique(dat$pop))){
    temp=cbind(sid_list[[i]],normalize_logscale[[i]]);      normalize_gm=rbind(normalize_gm,temp)}
    normalize_gm=data.frame(normalize_gm); names(normalize_gm)=c("Sample.ID","normal")
      return(normalize_gm)

}

rel_ancom<-function(OTUdat, Vardat, main.var="country", pr=NULL, pii=NULL,ref_name=NULL, paired=FALSE, test=TRUE){
  ########## STEP A1 ################################################################
  A1=struc_zero(OTUdat = OTUdat, Vardat = Vardat, main.var = main.var, p=pr)
  ########## STEP A2 ################################################################
  A2=reset_dat_missing(dat=A1[[1]],pii=pii,ref_name=ref_name,zeros=A1[[2]])
  adjusted_data=A2[[1]]; number_A2_adj_microbes=A2[[2]]; names_A2_adj_microbes=A2[[3]]
  ########## STEP A3 ################################################################
  adjusted_data[adjusted_data==0]=1
  ########## STEP B ################################################################
  ########## filter struc zero################################################################
  struc_zero_elim_ind= which(colSums(A1[[2]])>1);   sid=adjusted_data[,1];  
  normalu=adjusted_data[,3];   popu=adjusted_data[,2];  dat=adjusted_data[,-c(1,2,3)]
  ###Normalize and test individually 
  if(length(struc_zero_elim_ind)>0) redat=log(dat[,-struc_zero_elim_ind])-normalu else redat=log(dat)-normalu
  gauss_data=data.frame(Sample.ID=sid, pop=popu, redat)
  #############subset only the microbes#######################################
  gauss_microbiome=gauss_data %>% select(-c(Sample.ID,pop))
  cond <- apply(gauss_microbiome, MARGIN=2, FUN=function(x) if(any(is.na(x))) return(FALSE) else return(TRUE))
  gauss_microbiome <- gauss_microbiome[, which(cond==TRUE)]
  ########################analysis ##########################################
  d=dim(gauss_microbiome)[2]; pval=matrix(0,d,1); sig_ind=matrix(0,d,1)
  if(test==TRUE){
  for (j in 1:d){
	covariate=gauss_data$pop; response=gauss_microbiome[,j]
	regress=data.frame(response,covariate);  regress2=na.omit(regress)
	if(paired==FALSE){
		model=lm(regress2$response~ factor(regress2$covariate));  anv=anova(model);  pval[j]=anv$`Pr(>F)`[[1]]
		#pval[j]=t.test(regress2$response~ factor(regress2$covariate))$p.value
	} else {
		pval[j]=t.test(regress2$response~ factor(regress2$covariate), paired=TRUE)$p.value
	}
  }
  }
  padj=p.adjust(pval,"BH")
  #ind=which(padj<0.05)
  #detected_microbes=names(gauss_microbiome)[ind]
  df <- data.frame(feature=names(gauss_microbiome), ANCOMII.pval=pval, ANCOMII.adjpval.BH=padj, stringsAsFactors=FALSE)
  myord <- order(df$ANCOMII.pval, decreasing=FALSE)
  df <- df[myord, ]
  structural_zeros= data.frame(Microbe=names(adjusted_data[-c(1,2,3)]),t(1-A1[[2]]), stringsAsFactors=FALSE)
  return(list(df=df, structural_zeros=structural_zeros, gauss_data=gauss_data))

}






