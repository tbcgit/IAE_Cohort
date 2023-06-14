#!/usr/bin/env Rscript

suppressWarnings(suppressMessages(library(dplyr)))

#################################################
# Collapse rows in a data frame
#################################################
parameters <- commandArgs(trailingOnly = TRUE)

INPUTFILE  <- as.character(parameters[1])
GROUPCOL   <- as.character(parameters[2])
COLUMNS    <- as.character(parameters[3])
FUNCTION   <- as.character(parameters[4])
OUTPUTFILE <- as.character(parameters[5])

df <- read.table(file=INPUTFILE, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", 
                 check.names=FALSE, comment.char="", colClasses="character")

if(COLUMNS!="")
{
  columns <- c(GROUPCOL,
               sapply(strsplit(COLUMNS, split=",", fixed=TRUE)[[1]], 
                      trimws, USE.NAMES=FALSE))
  df <- df[,columns]
}                      

for(x in setdiff(colnames(df),GROUPCOL))
{
  df[,x] <- as.numeric(df[,x])
}

df <- df[order(df[,GROUPCOL]),]

g <- group_by(df, !!rlang::sym(GROUPCOL))
g <- summarise_all(g, rlang::as_function(FUNCTION))
df <- as.data.frame(g, stringsAsFactors=FALSE)

write.table(df, file=OUTPUTFILE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

