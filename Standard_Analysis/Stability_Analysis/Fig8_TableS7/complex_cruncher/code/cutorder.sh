#!/bin/sh

# ---------------------------------------------------
# Cut columns in a table 
# with the specified order of the columns
# ---------------------------------------------------

TABLE=$1
COLS=$2


for COL in $(echo $COLS | tr "," "\n")
do

  cut -f $COL $TABLE > aux_col.2.tsv
  
  if [ -f aux_tb.1.tsv ]
  then
     paste aux_tb.1.tsv aux_col.2.tsv > aux_tb.tsv
  else
     cp aux_col.2.tsv aux_tb.tsv
  fi

  cp aux_tb.tsv aux_tb.1.tsv

done


rm aux_tb.1.tsv
rm aux_col.2.tsv

cat aux_tb.tsv
rm aux_tb.tsv

