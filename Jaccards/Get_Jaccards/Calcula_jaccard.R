# 
wd <- '/media/benja/Extreme SSD/Tesi/DadesSusana/Treball_AttenuationBuffering/JaccardFiltBySamuel/IAE_16S'
setwd(wd)
dataRaw <- read.csv("TaxSilva_GEN_transformed.tsv", header=TRUE, sep='\t')
rownames(dataRaw) <- dataRaw$Genus
class(dataRaw)
colnames(dataRaw)
rownames(dataRaw)
sort(as.vector(colnames(dataRaw)))

library(dplyr)
dataOK <- select(dataRaw,-Genus)
Tdata <-t (dataOK)

#library(vegan)
#jaccard_table_distance <- vegdist(Tdata, method='jaccard', diag = TRUE, binary = FALSE, upper=TRUE)

library(jacpop)
jaccard_similarity <- generate_pw_jaccard(Tdata, pop.label = NULL, n.pcs = 10, plot_it = TRUE)
jaccard_similarity_df <- as.data.frame(jaccard_similarity$Jac)
colnames(jaccard_similarity_df) <- colnames(dataOK)
rownames(jaccard_similarity_df) <- colnames(dataOK)

write.csv(as.matrix(jaccard_similarity_df), file='jaccard_table_GENjacpop_fromNoNA.csv')


