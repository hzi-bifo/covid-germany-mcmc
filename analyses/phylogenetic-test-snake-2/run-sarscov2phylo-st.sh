#!/bin/bash

set -o xtrace
set -e
set -u # or set -o nounset

DIR=$1
cd $DIR
DATE_SEQ=$2
DATE_TREE=$3
X=$4


DIR=`pwd`
#cd $DIR/
cd results/fasttree/$DATE_TREE-$X/ ;
bash $DIR/sarscov2phylo/global_tree_gisaid.sh -i ./gisaid-$DATE_SEQ-$X-seq0.fa -o gisaid-$DATE_SEQ-ft -t 10
#conda deactivate sarscovphylo
cd $DIR
cp results/fasttree/$DATE_TREE-$X/ft_SH.tree results/gisaid-$DATE_TREE-$X-1.tree

