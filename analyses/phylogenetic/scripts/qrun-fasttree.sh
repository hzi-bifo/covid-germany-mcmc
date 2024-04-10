#!/bin/bash

# qsub -cwd -l h_vmem=10G,mem_free=10G,s_vmem=10G -pe smp 30 -o log -e log qrun.sh DATE_SEQ

conda activate sarscov2phylo

set -o xtrace
set -e
set -u # or set -o nounset

DATE_SEQ=$1
DATE_TREE=$2
X=$3
DIR=analyses/phylogenetic
DIR="$PWD"
cd $DIR/

cd results/fasttree/$DATE_TREE-$X/

bash $DIR/sarscov2phylo/global_tree_gisaid.sh -i ./gisaid-$DATE_SEQ-$X-seq0.fa -o gisaid-$DATE_SEQ-ft -t 10

cd $DIR
#cd ../../../
cp results/fasttree/$DATE_TREE-$X/ft_SH.tree results/gisaid-$DATE_TREE-$X-1.tree

