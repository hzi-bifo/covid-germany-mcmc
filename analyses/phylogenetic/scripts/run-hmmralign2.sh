#!/bin/bash

__conda_setup="$('/home/hforoughmand/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
#/home/hforoughmand/miniconda3/etc/profile.d/conda.sh
export PATH="/home/hforoughmand/miniconda3/bin:$PATH"

conda activate covid-uk

cd /net/viral_genomics/covid-lineage/huge-lineage-dynamics/analyses/phylogenetic

set -e
set -u # or set -o nounset
set -o xtrace

DATE_TREE=$1
X=$2
INDEX=$3

hmmalign2 --outformat A2M -o results/beast/unsampled/seq-split/$X-aln/$X-DTA-$DATE_TREE.hmm.$INDEX.aln results/beast/unsampled/$X-DTA-$DATE_TREE.0.hmm results/beast/unsampled/seq-split/$X/gisaid-$DATE_TREE-$X-unseq0.fa.$INDEX


echo "done"
