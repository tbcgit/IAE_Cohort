############
##

DATA=pathway.kegg.tpm.tsv
META=metaMTG.tsv
TAG=pathway.kegg.MTG

############ Sin agregar

./doTestDESeq2.R $DATA $META "group" "A,E,I" $TAG TRUE 0.6 0 DESCRIPTION FALSE TRUE "parametric"

./doTestDESeq2.R $DATA $META "group.time" "A.T1,A.T4,A.T5,A.T6,A.T7,A.T8" $TAG.timeEvolution.A TRUE 0.6 0 DESCRIPTION FALSE TRUE "parametric"
./doTestDESeq2.R $DATA $META "group.time" "E.T1,E.T4,E.T5,E.T6,E.T7,E.T8" $TAG.timeEvolution.E TRUE 0.6 0 DESCRIPTION FALSE TRUE "parametric"
./doTestDESeq2.R $DATA $META "group.time" "I.T1,I.T4,I.T5,I.T6,I.T7,I.T8" $TAG.timeEvolution.I TRUE 0.6 0 DESCRIPTION FALSE TRUE "parametric"

############ Agregando

./aggregateByFactor.R $DATA $META "donor" $TAG DESCRIPTION

./doTestDESeq2.R mean.counts.$TAG.tsv metadata.$TAG.tsv "group" "A,E,I" $TAG TRUE 0.6 0 DESCRIPTION FALSE TRUE "parametric"

