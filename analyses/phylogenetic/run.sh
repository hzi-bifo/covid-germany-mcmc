#DATE_TREE=$1
#DATE_METADATA=$2
#STATE=Germany

DATE_TREE=20210428
DATE_METADATA=20210503
DATE_SEQ=20210509
STATE=Germany

#V1: download tree from GISAID
#V2: build the tree
#V3: build tree as pangolin tree: 
# DATE_TREE=20210520
# DATE_METADATA=20210520
# DATE_SEQ=20210520
# ./scripts/tree-build-as-pangolin --out results/pangolin-$DATE_TREE-r0.tree --metadata ../../data/phylogenetic/gisaid-$DATE_METADATA-metadata.tsv --seqs ../../data/phylogenetic/gisaid-$DATE_SEQ.fa  -l Germany nonGermany
# ./scripts/phylogeo_sankoff_general --in results/pangolin-$DATE_TREE-r0.tree --out results/gisaid-$DATE_TREE-$STATE.tree --metadata ../../data/phylogenetic/gisaid-$DATE_METADATA-metadata.tsv --location_label $STATE non$STATE --cond 2 "==" Germany 



./scripts/phylogeo_sankoff_general --in ../../data/phylogenetic/gisaid-$DATE_TREE.tree --out results/gisaid-$DATE_TREE-$STATE.tree --metadata ../../data/phylogenetic/gisaid-$DATE_METADATA-metadata.tsv --location_label $STATE non$STATE --cond 2 "==" Germany --ilabel true
./scripts/contract_short_branch --in results/gisaid-$DATE_TREE-$STATE.tree --out results/gisaid-$DATE_TREE-cont.tree --short 5e-6 --location_label $STATE non$STATE --print-annotation true --print-internal-node-label true 

#V1:
./scripts/sampling --in results/gisaid-$DATE_TREE-cont.tree --out results/gisaid-$DATE_TREE-sampled.tree --short 5e-4 --shortstate 5e-6 -l $STATE non$STATE --samples-out results/gisaid-$DATE_TREE-samples.txt 
#V2:
./scripts/sample-evenly --in results/gisaid-$DATE_TREE-cont.tree --out results/gisaid-$DATE_TREE-sampled.tree~ -l $STATE non$STATE --samples-out results/gisaid-$DATE_TREE-samples.txt~ --metadata ../../data/phylogenetic/gisaid-$DATE_METADATA-metadata.tsv --bucket-size 100 800 --special 5 inode0 inode64 
# SUB_TREES="B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.1.214 B.1.177 B.1.160 B.1.221 B.1.36 B.1.258 B.1.351 P A C"
# ./scripts/sample-evenly --in results/gisaid-$DATE_TREE-cont.tree --out results/gisaid-$DATE_TREE-sampled.tree~ -l $STATE non$STATE --samples-out results/gisaid-$DATE_TREE-samples.txt~ --metadata ../../data/phylogenetic/gisaid-$DATE_METADATA-metadata.tsv --bucket-size 100 800 --special 5 $SUB_TREES --keep 1000 $SUB_TREES 
# cp results/gisaid-$DATE_TREE-sampled.tree~ results/gisaid-$DATE_TREE-sampled.tree
# cp results/gisaid-$DATE_TREE-samples.txt~ results/gisaid-$DATE_TREE-samples.txt

./scripts/tree-stat --in results/gisaid-$DATE_TREE-cont.tree -l $STATE non$STATE --large 1000 --depth 10
./scripts/tree-stat --in results/gisaid-$DATE_TREE-sampled.tree -l $STATE non$STATE --large 3000 --depth 30

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/"
./scripts/extract-metadata.R results/gisaid-$DATE_TREE-samples.txt ../../data/phylogenetic/gisaid-$DATE_METADATA-metadata.tsv results/gisaid-$DATE_TREE-metadata-sampled.tsv
# cut -f1 results/gisaid-$DATE_TREE-metadata-sampled.tsv | tail -n +2 | sed 's/^hCoV-19\///' > results/gisaid-$DATE_TREE-sample-names.txt 
# ./scripts/alignment-filter ../../data/phylogenetic/gisaid-$DATE_SEQ.fa.xz results/gisaid-$DATE_TREE-sample-names.txt > ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled0.fa 
# ./scripts/rename-alignment-ids results/gisaid-$DATE_TREE-metadata-sampled.tsv < ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled0.fa > ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled.fa
./scripts/alignment-filter ../../data/phylogenetic/gisaid-$DATE_SEQ.fa.xz results/gisaid-$DATE_TREE-samples.txt > ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled.fa 


SUB_TREES="inode18683 inode99244 inode239376 inode313849 inode199848"
./scripts/partition-by-name --in results/gisaid-$DATE_TREE-sampled.tree --par $SUB_TREES --samples "results/gisaid-$DATE_TREE-?-samples.txt" --trees "results/gisaid-$DATE_TREE-?-sampled.tree" -l $STATE non$STATE --print-annotation false --print-internal-node-label false
#Create sampled tree from the original tree (not the contracted one)
#./scripts/sampling-by-name --in results/gisaid-$DATE_TREE-$STATE.tree --out results/gisaid-$DATE_TREE-noncontract-sampled.tree -l $STATE non$STATE --samples-in results/gisaid-$DATE_TREE-samples.txt  (did not work, why?)
#./scripts/partition-by-name --in results/gisaid-$DATE_TREE-noncontract-sampled.tree --par inode18683 inode99244 inode239376 inode313849 inode199848 --samples "results/gisaid-$DATE_TREE-?-noncontract-samples.txt" --trees "results/gisaid-$DATE_TREE-?-noncontract-sampled.tree" -l $STATE non$STATE
./scripts/sample-and-partition-by-name --in results/gisaid-$DATE_TREE-$STATE.tree --samples "results/gisaid-$DATE_TREE-?-noncontract-samples.txt" --trees "results/gisaid-$DATE_TREE-?-noncontract-sampled.tree" -l $STATE non$STATE --samples-in "results/gisaid-$DATE_TREE-?-samples.txt" --par $SUB_TREES --print-annotation false --print-internal-node-label false

for X in $SUB_TREES; do 
	# ./scripts/fill-template.py --template data/X.fixedRootPrior.skygrid-template.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples.txt --tree results/gisaid-$DATE_TREE-$X-sampled.tree --tree_init results/gisaid-$DATE_TREE-$X-sampled.tree --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.xml --loc Germany --date $DATE_TREE 
	./scripts/fill-template.py --template data/X.fixedRootPrior.skygrid-template.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples.txt --tree results/gisaid-$DATE_TREE-$X-sampled.tree --tree_init results/gisaid-$DATE_TREE-$X-noncontract-sampled.tree --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.xml --loc Germany --date $DATE_TREE 
done

# run at grid.bifo
# DATE_TREE=20210520
# SUB_TREES="B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.1.214 B.1.177 B.1.160 B.1.221 B.1.36 B.1.258 B.1.351 P A C"
DATE_TREE=20210428
cd analyses/phylogenetic/
for i in {2..2}; do 
for X in $SUB_TREES; do 
qsub -cwd -l h_vmem=64G,mem_free=20G,s_vmem=20G -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast.sh run/$X-$i $DATE_TREE $X
done
done

for i in {1..1}; do 
for X in $SUB_TREES; do 
loganalyser results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-20210428.log
done
done

for i in {1..1}; do 
for X in $SUB_TREES; do 
python scripts/resample.py --burnin 15000000 --rate 100000 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-20210428.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-20210428-sampled.trees
python scripts/resample.py --burnin 15000000 --count 500 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-20210428.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-20210428-sub500.trees
done
done

#second rount of MCMC, DTA

for X in $SUB_TREES; do 
	./scripts/fill-template.py --template data/X-DTA-template.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples.txt --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X-DTA-$DATE_TREE.xml --loc Germany --date $DATE_TREE 
done

# run at grid.bifo
DATE_TREE=20210428
cd analyses/phylogenetic/
for i in {1..1}; do 
for X in $SUB_TREES; do 
qsub -cwd -l h_vmem=64G,mem_free=20G,s_vmem=20G -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA.sh run/$X-$i $DATE_TREE $X
done
done

for i in {1..1}; do 
for X in $SUB_TREES; do 
python scripts/resample.py --burnin 500000 --rate 4500 --tree results/beast/run/$X-$i/$X-DTA-20210428.trees --out results/beast/run/$X-$i/$X-DTA-20210428.sub4500.trees
done
done

for i in {1..1}; do 
for X in $SUB_TREES; do 
treeannotator -type mcc results/beast/run/$X-$i/$X-DTA-20210428.sub4500.trees results/beast/run/$X-$i/$X-DTA-20210428.MCC.trees
done
done

