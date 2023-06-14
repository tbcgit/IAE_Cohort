########
## Ver si hay diferencia para xW_V_stan entre los grupos definidos por Group, Gender, Gene, Mutation

./doTestOneVar.R xW_V_stan.tsv meta.tsv "group" "I,A,E" xW_V_stan_group FALSE 0 0 ID FALSE FALSE
./doTestOneVar.R xW_beta_stan.tsv meta.tsv "group" "I,A,E" xW_beta_stan_group FALSE 0 0 ID FALSE FALSE

./getPositiveMatrix.R xW_V_beta_stan.tsv meta.tsv xW_V_beta_stan_pos.tsv ID
./getEuclidianDist.R xW_V_beta_stan_pos.tsv meta.tsv mydist.csv
./doAdonisFromDistances.R mydist.csv meta.tsv "group" "I,A,E" adonis.V_beta TRUE FALSE 0

./getEuclidianDist.R RSI.matrix.tsv meta.tsv mydist.csv
./doPCOA.R mydist.csv meta.tsv "group" "I,A,E" pcoa.group TRUE FALSE 0.15
./doTest.R RSI.matrix.tsv meta.tsv "group" "I,A,E" group FALSE 0 0 ID FALSE FALSE

