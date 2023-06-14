
PROJECTDIR=/scratch/jpons/projects/susanaruiz/projectTimes
WORKDIR=$PROJECTDIR/jpons/20220809_kegg_complex_cruncher
ARTACHODIR=/scratch/aartacho/8-junio-2022-cohorte-IAE
SCRIPTSDIR=/scratch/jpons/scripts
PYTHONDIR=/scratch/software/python/bin


for METATYPE in "MTG" "MTT"
do
  echo $METATYPE

  mkdir $WORKDIR/$METATYPE 
  cd    $WORKDIR/$METATYPE
  
  # Matriz TPM de KEGGs filtrada
  # Si no filtramos filas, hay demasiados KEGGs y ComplexCruncher demanda demasiada RAM
  cat $ARTACHODIR/filtered.freqs.kegg.tpm.$METATYPE.tsv | sed "s/_$METATYPE//g" > KO.tpm.$METATYPE.filtered.tsv

  
  # Dividir por individuo
  mkdir $WORKDIR/$METATYPE/individuals
  cd    $WORKDIR/$METATYPE/individuals
  
  INDIVIDUALS=$(head -n 1 $WORKDIR/$METATYPE/KO.tpm.$METATYPE.filtered.tsv | cut -f1 --complement | sed "s/T[1-9]//g" | tr "\t" "\n" | sort | uniq)
  
  for IND in $INDIVIDUALS
  do
    echo $IND

    # Columnas correspondientes al individuo, ordenadas por tiempo
    INDCOLS=$(head -n 1 $WORKDIR/$METATYPE/KO.tpm.$METATYPE.filtered.tsv | tr "\t" "\n" | grep -n $IND"T[1-9]" | tr ":" "\t" | sort -k2 | cut -f1 | tr "\n" "," | rev | sed "s/,//" | rev)

    # Adultos como grupo de referencia sanos
    TBNAME=$IND
    AGE=$(echo $IND | cut -c 1)
    if [ $AGE == "A" ]; then TBNAME="h_$IND"; fi

    $SCRIPTSDIR/sh/cutorder.sh $WORKDIR/$METATYPE/KO.tpm.$METATYPE.filtered.tsv "1,$INDCOLS" | head -n 1 | sed "s/$IND//g"  > $TBNAME.tsv
    $SCRIPTSDIR/sh/cutorder.sh $WORKDIR/$METATYPE/KO.tpm.$METATYPE.filtered.tsv "1,$INDCOLS" | tail -n +2                  >> $TBNAME.tsv
  done


  # Generar Excel
  $PYTHONDIR/python $SCRIPTSDIR/python/tsvs2xlsx.py -i *.tsv -o $WORKDIR/$METATYPE/KO.tpm.$METATYPE.filtered.xlsx


  # TODO: Ejecutar Complex Cruncher: Máquina virtual


done



