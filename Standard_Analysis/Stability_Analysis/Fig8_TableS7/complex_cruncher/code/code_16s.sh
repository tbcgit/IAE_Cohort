PROJECTDIR=/scratch/jpons/projects/susanaruiz/projectTimes
WORKDIR=$PROJECTDIR/jpons/20220729_16s_complex_cruncher
SAMUELDIR=$PROJECTDIR/samuel/16S
SCRIPTSDIR=/scratch/jpons/scripts

cd $WORKDIR

# Resultados dada2
cp $SAMUELDIR/ASVs_counts_DEF.tsv .
cp $SAMUELDIR/ASVs_taxonomy_SILVA_v138_both_thr50_DEF_last.tsv .


# Agrupar por género
cut -f8 ASVs_taxonomy_SILVA_v138_both_thr50_DEF_last.tsv > a.txt
cut -f1 --complement ASVs_counts_DEF.tsv > b.tsv
paste a.txt b.tsv > counts.asv.tsv
rm a.txt b.tsv

$SCRIPTSDIR/r/summarise.table.R counts.asv.tsv "last" "" "sum" counts.genus.tsv
sed -i "s/^last/GENUS/" counts.genus.tsv
sed -i "s/16S_//g"      counts.genus.tsv


# Dividir por individuo
mkdir $WORKDIR/individuals
cd    $WORKDIR/individuals

INDIVIDUALS=$(head -n 1 $WORKDIR/counts.genus.tsv | cut -f1 --complement | sed "s/T[1-9]//g" | tr "\t" "\n" | sort | uniq)

for IND in $INDIVIDUALS
do
  echo $IND
 
  INDCOLS=$(head -n 1 $WORKDIR/counts.genus.tsv | tr "\t" "\n" | grep -n $IND"T[1-9]" | cut -f1 -d ":" | tr "\n" "," | rev | sed "s/,//" | rev)

  cut -f 1,$INDCOLS $WORKDIR/counts.genus.tsv | head -n 1 | sed "s/$IND//g"  > $IND.tsv
  cut -f 1,$INDCOLS $WORKDIR/counts.genus.tsv | tail -n +2                  >> $IND.tsv
done


# Generar Excel
cd $WORKDIR/individuals
$SCRIPTSDIR/python/tsvs2xlsx.py -i *.tsv -o $WORKDIR/counts.genus.individuals.xlsx

# TODO: Indicar pestañas healthy (3 excels)

# TODO: Ejecutar Complex Cruncher: Máquina virtual


