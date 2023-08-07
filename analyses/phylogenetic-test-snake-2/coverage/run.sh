conda activate r4-base

#Metadata from Hadi (containing a column with the importation lineage) was reformatted using metadata_reformatting.R
Rscript scripts/metadata_reformatting.R data/

#2.
#SDplots (https://github.com/hzi-bifo/SDplots_covlineages) were run using the reformatted metadata and the months file in data as input:

cd SDplots_covlineages/
bash SDPlots_lineages.sh ../data/clusterSamples_DTA_MCC_0.5_reformatted.tsv ../data/months.txt ../output 0.1
cd -

#3.
#Substitutions in significant lineages were investigated using lineage_substitutions.R
Rscript scripts/lineage_substitutions.R
