#!/bin/bash

#__conda_setup="$('/home/hforoughmand/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
#eval "$__conda_setup"
##/home/hforoughmand/miniconda3/etc/profile.d/conda.sh
#export PATH="/home/hforoughmand/miniconda3/bin:$PATH"
#
conda activate covid-uk

set -e
set -u # or set -o nounset
set -o xtrace

DIR=$1
DATE_TREE=$2
X=$3

#cd /net/viral_genomics/covid-lineage/huge-lineage-dynamics/analyses/phylogenetic/results/beast/
cd results/beast/
rm -rf $DIR
mkdir $DIR/
cp $X.fixedRootPrior.skygrid-$DATE_TREE.xml $DIR/
echo "executing beast"
cd $DIR/

#export PATH=/net/viral_genomics/covid-lineage/germany-lineage-dynamics/bin/beast/bin:$PATH
beast -save_every 1000000 -save_state $X.fixedRootPrior.skygrid-$DATE_TREE.state $X.fixedRootPrior.skygrid-$DATE_TREE.xml 

echo "done"
