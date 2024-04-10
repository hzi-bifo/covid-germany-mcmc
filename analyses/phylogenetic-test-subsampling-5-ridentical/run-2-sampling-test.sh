# In this file we do not add unsampled ones, 
# Testing differnet strategies of sampling with run-2.sh

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/"
#add beast into path: 
export PATH=$PATH_TO_BEAST/bin:$PATH
cd analyses/phylogenetic-test-subsampling-5-ridentical/
DATE_TREE=20210602
DATE_METADATA=20210602
DATE_SEQ=20210602
STATE=Germany
TWD=/tmp/
SUB_TREES="A B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.258 B.1.351 C"
FOLDER_MAIN=../phylogenetic/

OLD=../phylogenetic-test-subsampling-5
ln -s ../phylogenetic-test-snake-2/scripts .
mkdir -p results/beast/run/all/
for X in $SUB_TREES; do
	cp $OLD/results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree
done
cp $OLD/results/gisaid-$DATE_TREE-metadata-unsampled.tsv results/gisaid-$DATE_TREE-metadata-unsampled.tsv
cp $OLD/results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv

cp -r $OLD/reports-ius-3/ .
mkdir data
cp -r $OLD/data/*.csv data/

mkdir -p results/beast/unsampled/
mkdir -p results/beast/run/lin-ius/
mkdir -p results/beast/run/lin-ius/out/

mkdir -p temp/lin-ius-3/
cp $OLD/temp/lin-ius-3/*-2.log temp/lin-ius-3/
for X in $SUB_TREES; do
cat temp/lin-ius-3/$X-2.log | sed -n '/query_min_dist_lineage_half/,$p' | sed '/lin assigned/q' | sed '1d' | sed '$d' 
done | cut -d' ' -f3 | sort | uniq > temp/lin-ius-3/query-ids.txt
#
for X in $SUB_TREES; do
cat temp/lin-ius-3/$X-2.log | sed -n '/query_min_dist_lineage_half/,$p' | sed '/lin assigned/q' | sed '1d' | sed '$d' 
done | cut -d' ' -f4 | sort | uniq > temp/lin-ius-3/sample-idx.txt

for X in $SUB_TREES; do
cat temp/lin-ius-3/$X-2.log | sed -n '/query_min_dist_lineage_half/,$p' | sed '/lin assigned/q' | sed '1d' | sed '$d'
done > temp/lin-ius-3/match-mash-info.txt

for X in $SUB_TREES; do
cat temp/lin-ius-3/$X-2.log | sed -n '/query_min_dist_lineage_half/,$p' | sed '/lin assigned/q' | sed '1d' | sed '$d' 
done | cut -d' ' -f3,4  > temp/lin-ius-3/match-ids.txt

./scripts/alignment-filter $OLD/data/mmsa-$DATE_SEQ-sampled-unsampled0.fa <( cat temp/lin-ius-3/query-ids.txt temp/lin-ius-3/sample-idx.txt ) > temp/lin-ius-3/query-sample-seqs.fa 

for X in $SUB_TREES; do
	echo scripts/filter-alignment-distance $X
	echo -e '#!/bin/bash'"\nscripts/filter-alignment-distance <( cat temp/lin-ius-3/query-ids.txt temp/lin-ius-3/sample-idx.txt ) <( zcat $OLD/results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) > temp/lin-ius-3/alignment-distance-injected-$X.txt" > temp-filter-alignment-distance-$X.sh 
done 

for X in $SUB_TREES; do
sbatch -J "fa-$X" --qos verylong -t 24:00:00 temp-filter-alignment-distance-$X.sh
done

for X in $SUB_TREES; do
cat temp/lin-ius-3/alignment-distance-injected-$X.txt
done > temp/lin-ius-3/alignment-distance-injected.txt
rm temp-filter-alignment-distance-*.sh

#python scripts/compare-sequences-mash.py temp/lin-ius-3/alignment-distance-injected.txt temp/lin-ius-3/match-mash-info.txt > temp/lin-ius-3/sample-injected-distance-mash.txt 2> temp/lin-ius-3/sample-injected-distance-mash.err
python scripts/replace-mash-alignment-distance-with-real.py 0.0 temp/lin-ius-3/alignment-distance-injected.txt temp/lin-ius-3/query-sample-seqs.fa > temp/lin-ius-3/alignment-distance-real.txt 

#for X in B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.351 A C B.1.617.2; do
for X in $SUB_TREES; do
	echo "Tree $X "
	cp results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus
	./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation false --single-child false --ilabel true --set-tip-location false
	./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE-rich.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation true --single-child false --ilabel true --set-tip-location false
	
	sed 's/Germany+nonGermany/nonGermany/g' results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus > results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus
	mkdir -p results/beast/run/lin-ius-3/out/
	scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( cat temp/lin-ius-3/alignment-distance-real.txt ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius-3/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius-3/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage false --dist-threshold 0.0 --log results/beast/run/lin-ius-3/$X-2.log --do-assign true
done

function lin_ius_clusters_combine {
LINIUS_FOLDER=$1
(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/$LINIUS_FOLDER/out/$X-clusters_DTA_MCC_NA.tsv
for X in $SUB_TREES; do
	tail -n +2 results/beast/run/$LINIUS_FOLDER/out/$X-clusters_DTA_MCC_NA.tsv 
done ) > results/beast/run/$LINIUS_FOLDER/clusters_DTA_MCC_NA.tsv

(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/$LINIUS_FOLDER/out/$X-clusterSamples_DTA_MCC_NA.tsv
for X in $SUB_TREES; do
	tail -n +2 results/beast/run/$LINIUS_FOLDER/out/$X-clusterSamples_DTA_MCC_NA.tsv
done ) > results/beast/run/$LINIUS_FOLDER/clusterSamples_DTA_MCC_NA.tsv

(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/$LINIUS_FOLDER/out/$X-clusters_DTA_MCC_0.5.tsv
for X in $SUB_TREES; do
	tail -n +2 results/beast/run/$LINIUS_FOLDER/out/$X-clusters_DTA_MCC_0.5.tsv 
done ) > results/beast/run/$LINIUS_FOLDER/clusters_DTA_MCC_0.5.tsv

(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/$LINIUS_FOLDER/out/$X-clusterSamples_DTA_MCC_0.5.tsv
for X in $SUB_TREES; do
	tail -n +2 results/beast/run/$LINIUS_FOLDER/out/$X-clusterSamples_DTA_MCC_0.5.tsv
done ) > results/beast/run/$LINIUS_FOLDER/clusterSamples_DTA_MCC_0.5.tsv

}

#lin_ius_clusters_combine lin-ius
#lin_ius_clusters_combine lin-ius-2
lin_ius_clusters_combine lin-ius-3
#lin_ius_clusters_combine lin-ius-no
#lin_ius_clusters_combine lin-ius-3-aux

#git add results/beast/run/lin-ius/cluster*
git add results/beast/run/lin-ius-3/cluster*
git add data/*.csv 
git add results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv


