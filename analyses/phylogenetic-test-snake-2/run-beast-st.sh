#!/bin/bash

set -e
set -u # or set -o nounset
set -o xtrace

SDIR=$1
DIR=$2
DATE_TREE=$3
X=$4

cd $SDIR
rm -rf $DIR
mkdir $DIR/
cp $X.fixedRootPrior.skygrid-$DATE_TREE.xml $DIR/
echo "executing beast"
cd $DIR/

beast-v1.10.5pre_thorney_0.1.1 -save_every 1000000 -save_state $X.fixedRootPrior.skygrid-$DATE_TREE.state $X.fixedRootPrior.skygrid-$DATE_TREE.xml 

echo "done"
