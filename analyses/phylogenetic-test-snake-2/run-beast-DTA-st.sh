#!/bin/bash

set -e
set -u # or set -o nounset
set -o xtrace

DIR=$1
DATE_TREE=$2
X=$3

cd analyses/phylogenetic-test-snake-2/results/beast/
#rm -rf $DIR
mkdir -p $DIR/
cp $X-DTA-$DATE_TREE.xml $DIR/
# I added this cp! is it the correct version?
cp $X.fixedRootPrior.skygrid-$DATE_TREE.trees $DIR/
echo "executing beast"
cd $DIR/

#add beast into path: 
export PATH=$PATH_TO_BEAST/bin:$PATH
beast-c8cc55d4fe -save_every 1000000 -overwrite -save_state $X-DTA-$DATE_TREE.state $X-DTA-$DATE_TREE.xml

echo "done"
