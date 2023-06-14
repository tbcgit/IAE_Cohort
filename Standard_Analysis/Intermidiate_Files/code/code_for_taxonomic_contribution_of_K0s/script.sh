##############
## Cuidado!! Coger los scripts de la carpeta de trabajo y no de downstream.analysis, porque han sido editados ad hoc para este trabajo ...

cd /scratch/aartacho/5-julio-2022-cohorte-IAE-contribucion-taxonomica-KOs/species

for file in `ls K*.DNA.species.contribution.tsv`
do
../splitTaxLevels_conNAs.R $file ../metaMTG.tsv "k_;p_;c_;o_;f_;g_;s_" "sum"
echo $file
done

for file in `ls K*.RNA.species.contribution.tsv`
do
../splitTaxLevels_conNAs.R $file ../metaMTT.tsv "k_;p_;c_;o_;f_;g_;s_" "sum"
echo $file
done

rm s_*

##############

for file in `ls *DNA*`
do
study=`echo $file | sed 's/.tsv//g'`
../aggregateByFactor.R $file ../metaMTG.tsv "group" $study "Taxon"
echo $file
done

for file in `ls K*DNA*`
do
study=`echo $file | sed 's/.tsv//g'`
../aggregateByFactor.R $file ../metaMTG.tsv "group" $study "ID"
echo $file
done

for file in `ls *RNA*`
do
study=`echo $file | sed 's/.tsv//g'`
../aggregateByFactor.R $file ../metaMTT.tsv "group" $study "Taxon"
echo $file
done

for file in `ls K*RNA*`
do
study=`echo $file | sed 's/.tsv//g'`
../aggregateByFactor.R $file ../metaMTT.tsv "group" $study "ID"
echo $file
done

##############

for file in `ls p_*DNA* o_*DNA* c_*DNA* f_*DNA* g_*DNA* K*DNA*`
do
study=`echo $file | sed 's/.tsv//g'`
../getNormTable.R $file ../metaMTG.tsv $study.pct.tsv
echo $file
done

for file in `ls p_*RNA* o_*RNA* c_*RNA* f_*RNA* g_*RNA* K*RNA*`
do
study=`echo $file | sed 's/.tsv//g'`
../getNormTable.R $file ../metaMTT.tsv $study.pct.tsv
echo $file
done

for file in `ls agregado*`
do
study=`echo $file | sed 's/.tsv//g'`
../getNormTable.R $file ../meta_agregado.tsv $study.pct.tsv
echo $file
done

##############

mv p_* ../phylum/
mv agregado.p_* ../phylum/

mv o_* ../order/
mv agregado.o_* ../order/

mv c_* ../class/
mv agregado.c_* ../class/

mv f_* ../family/
mv agregado.f_* ../family/

mv g_* ../genus/
mv agregado.g_* ../genus/



