############
##

RAWDATA=ASVs_counts.tsv
META=meta16S.tsv

./doDiversityColorByGroup.R $RAWDATA $META asv 100 "group" "I,A,E"

############

./doANCOMBC_MT2G.R $RAWDATA $META "group" "A,E,I" sinfilt FALSE 0 0 TAXONOMY "group" FALSE

./getPositiveMatrix.R ANCOMBC.normalized.$RAWDATA $META ANCOMBC.norm.pos.$RAWDATA DESCRIPTION

DATA=ANCOMBC.norm.pos.$RAWDATA

############ Sin agregar

./getBrayCurtisDist.R $DATA $META dist.csv
./doPCOA.R dist.csv $META "group" "A,E,I" pcoa.asv TRUE FALSE 0

./doTest.R $DATA $META "group" "A,E,I" asv TRUE 0.6 0 DESCRIPTION FALSE FALSE

./doTest.R $DATA $META "group.time" "A.T1,A.T2,A.T3,A.T4,A.T5,A.T6,A.T7,A.T8" asv.timeEvolution.A TRUE 0.6 0 DESCRIPTION FALSE FALSE
./doTest.R $DATA $META "group.time" "E.T1,E.T2,E.T3,E.T4,E.T5,E.T6,E.T7,E.T8" asv.timeEvolution.E TRUE 0.6 0 DESCRIPTION FALSE FALSE
./doTest.R $DATA $META "group.time" "I.T1,I.T2,I.T3,I.T4,I.T5,I.T6,I.T7,I.T8" asv.timeEvolution.I TRUE 0.6 0 DESCRIPTION FALSE FALSE

./doANCOMBC.R $RAWDATA $META "group" "I,A" asv.IA TRUE 0.6 0 TAXONOMY "group" FALSE
./doANCOMBC.R $RAWDATA $META "group" "I,E" asv.IE TRUE 0.6 0 TAXONOMY "group" FALSE
./doANCOMBC.R $RAWDATA $META "group" "A,E" asv.AE TRUE 0.6 0 TAXONOMY "group" FALSE

############ Agregando

./aggregateByFactorNormData.R $DATA $META "donor" asv DESCRIPTION

./doTest.R averaged.asv.tsv averaged.metadata.asv.tsv "group" "A,E,I" asv TRUE 0.6 0 DESCRIPTION FALSE FALSE

