#!/bin/bash -e

set -e
set -u # or set -o nounset
set -o xtrace


DATE_TREE=$1
STATE=$2
TREE_NAME=$3

echo scripts/generate-effectiveness.sh 

#treeannotator -type mcc results/beast/run/$X-$i/$X-DTA-$DATE_TREE.sub4500.trees results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree
treeannotator -type mcc results/beast/run/lin-ius-3e/sample-tree/$TREE_NAME.nexus results/beast/run/lin-ius-3e/mcc-tree/$TREE_NAME.mcc.nexus

scripts/lineage-importation-inject-unsampled --tree results/beast/run/lin-ius-3e/mcc-tree/$TREE_NAME.mcc.nexus --dist /dev/null -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius-3e/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$TREE_NAME-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/run/lin-ius-3e/out-tree/$TREE_NAME.nexus --out-clusters results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius-3e/out/$TREE_NAME-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3e/out/$TREE_NAME-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage false --dist-threshold 0.0 --log results/beast/run/lin-ius-3e/log/inj-$TREE_NAME-2.log --do-assign true
Rscript scripts/generate-effectiveness.R results/beast/run/lin-ius-3e/out/$TREE_NAME- results/beast/run/lin-ius-3e/out/$TREE_NAME-effectiveness.tsv

echo "done"
