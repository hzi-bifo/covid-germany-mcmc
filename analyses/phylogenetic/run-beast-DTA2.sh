#!/bin/bash

set -e
set -u # or set -o nounset
set -o xtrace

DIR=$1
DATE_TREE=$2
X=$3

cd results/beast/
#rm -rf $DIR
#mkdir $DIR/
cp $X-DTA2-$DATE_TREE.xml $DIR/
echo "executing beast"
cd $DIR/

beast -save_every 1000000 -overwrite -save_state $X-DTA2-$DATE_TREE.state $X-DTA2-$DATE_TREE.xml

echo "done"
