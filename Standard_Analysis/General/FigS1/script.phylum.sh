############
##

RAWDATA=p_counts.tsv
META=meta16S.tsv

./aggregateByFactorNormDataSum.R $RAWDATA $META "group" phylum Taxon
./getNormTable.R averaged.phylum.tsv averaged.metadata.phylum.tsv phylum.group.aggregated.pct.tsv 1 "ID" FALSE

