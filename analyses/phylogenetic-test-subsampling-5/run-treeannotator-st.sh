#!/bin/bash

__conda_setup="$('/home/hforoughmand/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
#/home/hforoughmand/miniconda3/etc/profile.d/conda.sh
export PATH="/home/hforoughmand/miniconda3/bin:$PATH"

conda activate covid-uk

set -e
set -u # or set -o nounset
set -o xtrace

DATE_TREE=$1
X=$2
i=$3


echo "Treeannotator for $DATE_TREE $X $i"

cd /net/viral_genomics/covid-lineage/huge-lineage-dynamics/analyses/phylogenetic-test-subsampling-5/
export PATH=/net/viral_genomics/covid-lineage/germany-lineage-dynamics/bin/beast/bin:$PATH
#treeannotator -type mcc results/beast/run/$X-$i/$X-DTA-$DATE_TREE.sub4500.trees results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree
treeannotator -type mcc  results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees  results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree

echo "done"
