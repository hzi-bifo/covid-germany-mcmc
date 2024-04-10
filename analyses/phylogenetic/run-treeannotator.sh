#!/bin/bash

set -e
set -u # or set -o nounset
set -o xtrace

DATE_TREE=$1
X=$2
i=$3


echo "Treeannotator for $DATE_TREE $X $i"

treeannotator -type mcc  results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees  results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree

echo "done"
