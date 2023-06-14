#############
##

./doBarplot.R agregado.g_K01667.DNAcontribution.tsv meta.tsv "samples" 0.1 "ID" "genus"
mv summary.taxa.barplot.all.samples.pdf barplot.K01667.DNA.pdf

./doBarplot.R agregado.g_K01667.RNAcontribution.tsv meta.tsv "samples" 0.1 "ID"
mv summary.taxa.barplot.all.samples.pdf barplot.K01667.RNA.pdf

./doBarplot.R agregado.g_K01696.DNAcontribution.tsv meta.tsv "samples" 0.1 "ID"
mv summary.taxa.barplot.all.samples.pdf barplot.K01696.DNA.pdf

./doBarplot.R agregado.g_K01696.RNAcontribution.tsv meta.tsv "samples" 0.1 "ID"
mv summary.taxa.barplot.all.samples.pdf barplot.K01696.RNA.pdf

#######

./doBarplot.R g_K01667.tsv meta1.tsv "samples" 0.1 "ID"
mv summary.taxa.barplot.all.samples.pdf barplot.K01667.pdf

./doBarplot.R g_K01696.tsv meta1.tsv "samples" 0.1 "ID"
mv summary.taxa.barplot.all.samples.pdf barplot.K01696.pdf

#######

./doBarplot.R g_counts.tsv meta2.tsv "samples" 0.1 "ID"
mv summary.taxa.barplot.all.samples.pdf barplot.all.pdf



