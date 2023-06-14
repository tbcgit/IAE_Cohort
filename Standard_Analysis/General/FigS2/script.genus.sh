############
##

RAWDATA=g_counts.tsv
META=meta16S.tsv

./aggregateByFactorNormDataSum.R $RAWDATA $META "group" genus Taxon
./getNormTable.R averaged.genus.tsv averaged.metadata.genus.tsv genus.group.aggregated.pct.tsv 1 "ID" FALSE
./doBarplot.R averaged.genus.tsv averaged.metadata.genus.tsv "group" 0.7 "ID"

############

./doANCOMBC_MT2G.R $RAWDATA $META "group" "A,E,I" sinfilt FALSE 0 0 Taxon "group" FALSE

./getPositiveMatrix.R ANCOMBC.normalized.$RAWDATA $META ANCOMBC.norm.pos.$RAWDATA ID

DATA=ANCOMBC.norm.pos.$RAWDATA

############ Sin agregar

./getBrayCurtisDist.R $DATA $META dist.csv
./doPCOA.R dist.csv $META "group" "A,E,I" pcoa.genus TRUE FALSE 0

./doTest.R $DATA $META "group" "A,E,I" genus TRUE 0.6 0 ID FALSE FALSE

./doTest.R $DATA $META "group.time" "A.T1,A.T2,A.T3,A.T4,A.T5,A.T6,A.T7,A.T8" genus.timeEvolution.A TRUE 0.6 0 ID FALSE FALSE
./doTest.R $DATA $META "group.time" "E.T1,E.T2,E.T3,E.T4,E.T5,E.T6,E.T7,E.T8" genus.timeEvolution.E TRUE 0.6 0 ID FALSE FALSE
./doTest.R $DATA $META "group.time" "I.T1,I.T2,I.T3,I.T4,I.T5,I.T6,I.T7,I.T8" genus.timeEvolution.I TRUE 0.6 0 ID FALSE FALSE

./doANCOMBC.R $RAWDATA $META "group" "I,A" genus.IA TRUE 0.6 0 Taxon "group" FALSE
./doANCOMBC.R $RAWDATA $META "group" "I,E" genus.IE TRUE 0.6 0 Taxon "group" FALSE
./doANCOMBC.R $RAWDATA $META "group" "A,E" genus.AE TRUE 0.6 0 Taxon "group" FALSE

############ Agregando

./aggregateByFactorNormData.R $DATA $META "donor" genus ID

./doTest.R averaged.genus.tsv averaged.metadata.genus.tsv "group" "A,E,I" genus TRUE 0.6 0 DESCRIPTION FALSE FALSE

