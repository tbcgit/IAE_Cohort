First you will need to clone the data and src directories of the Eng & Borenstein (2018) pipeline
(https://github.com/borenstein-lab/robustness).

You will also need to add to the "data" directory the Picrust2 result data files:
- KO_predicted.tsv.gz
- marker_predicted_and_nsti.tsv.gz
- out.tre

Then the script ParseCounts_V3.py is used to generate count tables for each sample based on the 16S count table.
You will need to save the results in a directory called "TSV".

After that, the script OrdresRobustness.py is used to run the src code from Eng & Borenstein (2018) for each sample.

The result of the execution is: 
- On the one hand, one file per sample containing only the global Attenuation and Buffering data.
- On the other hand, another file per sample containing Attenuation and Buffering for certain superpathways and pathways.

In all cases this files are converted into a single table each using bash (to generate global values, values per superpathway and values per pathway):
----
grep 16S_ *robustness_factors | cut -d ":" -f 2 > TableRobustnessFactors_R20220220.tsv
sed -i '1i Sample\t\Attenuation\Buffering' TableRobustnessFactors_R20220220.tsv

grep 16S_ *_specific_robustness_factors | cut -d ":" -f 2 > TableSpecifRobustnessFactors_bySuperpathway.tsv
sed -i '1i original_sample\t\funct\t\attenuation\buffering' TableSpecifRobustnessFactors_bySuperpathway.tsv

grep ko *_specific_robustness_factors | cut -d ":" -f 2 >  TableAttBuff_by_ko_2paper.tsv
sed -i '1i Sample\t\Pathway\t\Attenuation\Buffering' TableAttBuff_by_ko_2paper.tsv

----

The final result tables are:
-TableSpecifRobustnessFactors_bySuperpathway.tsv [Values per superpathway]
-TableRobustnessFactors_R20220220.tsv [Global values]
-TableAttBuff_by_ko_2paper.xlsx (after converting "TableAttBuff_by_ko_2paper.tsv" to excel file) [Values per pathway]

That are then used to generate the final graphics and analysis:
- Figure 9
- Supplementary Figure 5
- Supplementary Figure 6
- Supplementary Table 8
- Supplementary Table 9
