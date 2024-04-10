# Testing differnet strategies of sampling with run-2.sh

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/"
#add beast into path: 
export PATH=$PATH_TO_BEAST/bin:$PATH
cd analyses/phylogenetic-test-snake-2/
DATE_TREE=20210602
DATE_METADATA=20210602
DATE_SEQ=20210602
STATE=Germany
TWD=/tmp/
SUB_TREES="A B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.258 B.1.351 C"
FOLDER_MAIN=../phylogenetic/

##cleaning (is done via run-2.sh, we skip it here)
#grep -A1000 -F 'hCoV-19/Wuhan/WH04/2020' ../../../data/data/gisaid-$DATE_SEQ-raw.fa > ../../../data/data/gisaid-$DATE_SEQ-Wu1.fa
##edit then ..., label: hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05
#DIR=`pwd`
#mkdir $TWD/gisaid-$DATE_SEQ/
#cd $TWD/gisaid-$DATE_SEQ/
#ln -s $DIR/../../../data/data/gisaid-$DATE_SEQ-raw.fa .
#bash $DIR/sarscov2phylo/clean_gisaid.sh -i gisaid-$DATE_SEQ-raw.fa -o gisaid-$DATE_SEQ.fa -t 30 
#mv gisaid-$DATE_SEQ.fa $DIR/../../../data/data/
#
#cd $DIR
#
### SUBSAMPLING: 
##V3 (CURRENT): build tree as pangolin tree: 
#./scripts/tree-build-as-pangolin --out results/pangolin-$DATE_TREE-r0.tree --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv --seqs ../../../data/data/gisaid-$DATE_SEQ.fa  -l Germany nonGermany
## we use following for filtering non-good samaple (e.g. sample without complete date)
#./scripts/phylogeo_sankoff_general --in results/pangolin-$DATE_TREE-r0.tree --out results/gisaid-$DATE_TREE-$STATE.tree --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv --location_label $STATE non$STATE --cond 2 "==" Germany 
##./scripts/contract_short_branch --in results/gisaid-$DATE_TREE-$STATE.tree --out results/gisaid-$DATE_TREE-cont.tree --short 5e-6 --location_label $STATE non$STATE --print-annotation true --print-internal-node-label true 

mkdir results/ results/beast
mkdir data/
ln -s ../$FOLDER_MAIN/results/gisaid-$DATE_TREE-$STATE.tree results/

Rscript scripts/filter-outliers.R '../../../data/lineage-list.csv' ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv results/gisaid-$DATE_METADATA-metadata.tsv

#TEST: ./scripts/sample-evenly --in results/gisaid-$DATE_TREE-$STATE.tree --out $TMP/results_gisaid-$DATE_TREE-sampled.tree~ -l $STATE non$STATE --samples-out $TMP/results_gisaid-$DATE_TREE-samples.txt~ --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv --bucket-size 100 15 --special 5 $SUB_TREES --keep 1000 10000 $SUB_TREES --seed 7
./scripts/sample-evenly --in results/gisaid-$DATE_TREE-$STATE.tree --out results/gisaid-$DATE_TREE-sampled.tree~ -l $STATE non$STATE --samples-out results/gisaid-$DATE_TREE-samples.txt~ --metadata results/gisaid-$DATE_METADATA-metadata.tsv --bucket-size 100 25 --special 5 $SUB_TREES --keep 0 10000 $SUB_TREES --unbias-method threshold 
# ./scripts/sample-evenly --in results/gisaid-$DATE_TREE-$STATE.tree --out temp/results/gisaid-$DATE_TREE-sampled.tree~ -l $STATE non$STATE --samples-out temp/results/gisaid-$DATE_TREE-samples.txt~ --metadata results/gisaid-$DATE_METADATA-metadata.tsv --bucket-size 100 25 --special 5 $SUB_TREES --keep 0 10000 $SUB_TREES --unbias-method threshold # TEST!
# ./scripts/sample-evenly --in results/gisaid-$DATE_TREE-$STATE.tree --out temp/results/gisaid-$DATE_TREE-sampled.tree~2 -l $STATE non$STATE --samples-out temp/results/gisaid-$DATE_TREE-samples.txt~2 --metadata results/gisaid-$DATE_METADATA-metadata.tsv --bucket-size 100 25 --special 5 $SUB_TREES --keep 0 10000 $SUB_TREES --unbias-method threshold # TEST!
cp results/gisaid-$DATE_TREE-sampled.tree~ results/gisaid-$DATE_TREE-sampled.tree
cp results/gisaid-$DATE_TREE-samples.txt~ results/gisaid-$DATE_TREE-samples.txt

# Trees are not there, why?
## We can check the statistics of the trees and subtrees with following commands. The SUB_TREES variable is calculated beased on following statistics.
#./scripts/tree-stat --in results/gisaid-$DATE_TREE.tree -l $STATE non$STATE --large 1000 --depth 10
#./scripts/tree-stat --in results/gisaid-$DATE_TREE-sampled.tree -l $STATE non$STATE --large 500 --depth 30

./scripts/extract-metadata.R results/gisaid-$DATE_TREE-samples.txt results/gisaid-$DATE_METADATA-metadata.tsv results/gisaid-$DATE_TREE-metadata-sampled.tsv

cut -f1 results/gisaid-$DATE_TREE-metadata-sampled.tsv | tail -n +2  > results/gisaid-$DATE_TREE-sample-names.txt 
# changed to local folders
./scripts/alignment-filter ../../../data/data/gisaid-$DATE_SEQ.fa results/gisaid-$DATE_TREE-sample-names.txt > data/gisaid-$DATE_SEQ-sampled0.fa 
grep '>' data/gisaid-$DATE_SEQ-sampled0.fa | sort | uniq -c | grep -v '^ *1 ' | sort -nr 
./scripts/rename-alignment-ids results/gisaid-$DATE_TREE-metadata-sampled.tsv < data/gisaid-$DATE_SEQ-sampled0.fa > data/gisaid-$DATE_SEQ-sampled.fa

# Extracting data of the SUB_TREES
./scripts/partition-by-name --in results/gisaid-$DATE_TREE-sampled.tree --par $SUB_TREES --samples "results/gisaid-$DATE_TREE-?-samples0.txt" --trees "results/gisaid-$DATE_TREE-?-sampled0.tree" -l $STATE non$STATE --print-annotation false --print-internal-node-label false
#statistics:
#for X in $SUB_TREES; do echo $X; 
#./scripts/phylogeo_sankoff_general --in results/gisaid-$DATE_TREE-$STATE.tree --out x.tree --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --location_label $STATE non$STATE --cond 2 "==" Germany 
#./scripts/tree-stat --in x.tree -l $STATE non$STATE --large 0 --depth 0; done
#./scripts/sample-and-partition-by-name --in results/gisaid-$DATE_TREE-$STATE.tree --samples "results/gisaid-$DATE_TREE-?-noncontract-samples.txt" --trees "results/gisaid-$DATE_TREE-?-noncontract-sampled.tree" -l $STATE non$STATE --samples-in "results/gisaid-$DATE_TREE-?-samples.txt" --par $SUB_TREES --print-annotation false --print-internal-node-label false

#create fasttree for each partition:
# Creating sequences for each sub-tree
for X in $SUB_TREES; do
	./scripts/alignment-filter data/gisaid-$DATE_SEQ-sampled.fa results/gisaid-$DATE_TREE-$X-samples0.txt > results/gisaid-$DATE_TREE-$X-seq0.fa
	cat ../../../data/data/gisaid-$DATE_SEQ-Wu1.fa >> results/gisaid-$DATE_TREE-$X-seq0.fa
	#cd results/fasttree/$DATE_TREE-$X/
	#bash $DIR/sarscov2phylo/global_tree_gisaid.sh -i gisaid-$DATE_TREE-$X-seq0.fa -o gisaid-$DATE_TREE-$X-tree -t 20
	#qsub -cwd -l h_vmem=10G,mem_free=10G,s_vmem=10G -pe smp 30 -o log -e log qrun.sh $DATE_SEQ $DATE_TREE $X
	#cd -
done


#run on grid
for X in $SUB_TREES; do
	DIR=`pwd` ;
	rm -rf results/fasttree/$DATE_TREE-$X/ ;
	mkdir -p results/fasttree/$DATE_TREE-$X/ ;
	ln -s $DIR/results/gisaid-$DATE_TREE-$X-seq0.fa results/fasttree/$DATE_TREE-$X/gisaid-$DATE_TREE-$X-seq0.fa ;
	#cp scripts/qrun-fasttree.sh results/fasttree/$DATE_TREE-$X/qrun.sh ;

	cat > results/fasttree/$DATE_TREE-$X/qrun.sh << EOF
#!/bin/bash

set -o xtrace
set -e
set -u # or set -o nounset

cd $DIR/
cd results/fasttree/$DATE_TREE-$X/ ;
bash $DIR/sarscov2phylo/global_tree_gisaid.sh -i ./gisaid-$DATE_SEQ-$X-seq0.fa -o gisaid-$DATE_SEQ-ft -t 21
#conda deactivate sarscovphylo
cd $DIR
rm -f results/gisaid-$DATE_TREE-$X-1.tree
cp results/fasttree/$DATE_TREE-$X/ft_FBP.tree results/gisaid-$DATE_TREE-$X-1.tree

EOF

	qsub -cwd -N sarscov2phylo-$X -M foroughmand@gmail.com -l h_vmem=30G,mem_free=30G,s_vmem=30G -pe smp 1 -o log -e log results/fasttree/$DATE_TREE-$X/qrun.sh $DATE_SEQ $DATE_TREE $X ;
	#JOB_ID=$RANDOM
	#bash results/fasttree/$DATE_TREE-$X/qrun.sh $DATE_SEQ $DATE_TREE $X 2>log/sarscov2phylo-$X.e$JOB_ID >log/sarscov2phylo-$X.o$JOB_ID
	

	#bash qrun.sh $DATE_SEQ $DATE_TREE $X
# qsub -cwd -l h_vmem=64G,mem_free=20G,s_vmem=20G -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA.sh run/$X-$i $DATE_TREE $X# qsub -cwd -l h_vmem=64G,mem_free=20G,s_vmem=20G -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA.sh run/$X-$i $DATE_TREE $X
# qsub -cwd -l h_vmem=64G,mem_free=20G,s_vmem=20G -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA.sh run/$X-$i $DATE_TREE $X
# -pe smp 10 
	#cd - ;
done

#Status
for X in $SUB_TREES; do 
GRID_ID=`ls log/sarscov2phylo-$X.e* | sed 's/.*e//' | sort | tail -n1`
LAST_LINE=`tail -n1 log/sarscov2phylo-$X.o$GRID_ID`
if [ "$LAST_LINE" != "done" ]; then
tail -n 3 log/sarscov2phylo-$X.o$GRID_ID log/sarscov2phylo-$X.e$GRID_ID
fi
done

declare -A SUB_TREE_CLOCK_RATE=(
		[B.1.1.7]=0.00075
		[B.1.1.317]=0.00075
		[B.1.1.214]=0.00075
		[B.1.160]=0.00075
		[B.1.221]=0.00075
		[B.1.36]=0.00075
		[B.1.351]=0.00075
		[P]=0.00075
		[A]=0.00075
		[C]=0.00075 
		[B.1.617.2]=0.00075
		[B.1.1.519]=0.00075
		[B.1.258]=0.00075
		[B.1.1.70]=0.00075
		[B.1.177]=0.00075
)
#remove case EPI_ISL_2333117 which has date 1905-06-26
#R data <- read.table('results/gisaid-20210602-metadata-sampled.tsv', sep='\t', quote='"', head = TRUE)
#  data[as.Date(data$Collection.date) <= "2019-10-22","Accession.ID"]
# OLD SAMPLES:  EPI_ISL_412860 EPI_ISL_402131 EPI_ISL_412976 EPI_ISL_412977 EPI_ISL_852605 EPI_ISL_852604 EPI_ISL_610156 EPI_ISL_1699443 EPI_ISL_2333072 EPI_ISL_2333073 EPI_ISL_2333074 EPI_ISL_2333075 EPI_ISL_2333082 EPI_ISL_2333083 EPI_ISL_2333084 EPI_ISL_2333086 EPI_ISL_2333087 EPI_ISL_2333088 EPI_ISL_2333089 EPI_ISL_2333090 EPI_ISL_2333091 EPI_ISL_2333092 EPI_ISL_2333093 EPI_ISL_2333094 EPI_ISL_2333095 EPI_ISL_2333096 EPI_ISL_2333097 EPI_ISL_2333098 EPI_ISL_2333099 EPI_ISL_2333100 EPI_ISL_2333101 EPI_ISL_2333102 EPI_ISL_2333103 EPI_ISL_2333104 EPI_ISL_2333105 EPI_ISL_2333106 EPI_ISL_2333107 EPI_ISL_2333108 EPI_ISL_2333109 EPI_ISL_2333110 EPI_ISL_2333111 EPI_ISL_2333112 EPI_ISL_2333113 EPI_ISL_2333114 EPI_ISL_2333115 EPI_ISL_2333116 EPI_ISL_2333117 EPI_ISL_2333178 EPI_ISL_2333182 EPI_ISL_2333183 EPI_ISL_2333184 EPI_ISL_2333185 EPI_ISL_2333186 EPI_ISL_2333187 EPI_ISL_2333190 EPI_ISL_2333191 EPI_ISL_2333192 EPI_ISL_2333193 EPI_ISL_2333194 EPI_ISL_2333195 EPI_ISL_2333196 EPI_ISL_2333202 EPI_ISL_2333203 EPI_ISL_2333214 EPI_ISL_2333215 EPI_ISL_2333216 EPI_ISL_2333217 EPI_ISL_2333516 EPI_ISL_2333517 EPI_ISL_2333518 EPI_ISL_2333519 EPI_ISL_2333520 EPI_ISL_2333521 EPI_ISL_2333522 EPI_ISL_2333523 EPI_ISL_2333524 EPI_ISL_2333525 EPI_ISL_2333526 EPI_ISL_2333527 EPI_ISL_2333528 EPI_ISL_2357883 

for X in $SUB_TREES; do
	#cp results/fasttree/$DATE_TREE-$X/ft_FBP.tree results/gisaid-$DATE_TREE-$X-1.tree ;
	./scripts/extract-metadata.R results/gisaid-$DATE_TREE-$X-samples0.txt results/gisaid-$DATE_TREE-metadata-sampled.tsv results/gisaid-$DATE_TREE-$X-metadata-sampled-1.tsv

	grep -v "EPI_ISL_412860\|EPI_ISL_402131\|EPI_ISL_412976\|EPI_ISL_412977\|EPI_ISL_852605\|EPI_ISL_852604\|EPI_ISL_610156\|EPI_ISL_1699443\|EPI_ISL_2333072\|EPI_ISL_2333073\|EPI_ISL_2333074\|EPI_ISL_2333075\|EPI_ISL_2333082\|EPI_ISL_2333083\|EPI_ISL_2333084\|EPI_ISL_2333086\|EPI_ISL_2333087\|EPI_ISL_2333088\|EPI_ISL_2333089\|EPI_ISL_2333090\|EPI_ISL_2333091\|EPI_ISL_2333092\|EPI_ISL_2333093\|EPI_ISL_2333094\|EPI_ISL_2333095\|EPI_ISL_2333096\|EPI_ISL_2333097\|EPI_ISL_2333098\|EPI_ISL_2333099\|EPI_ISL_2333100\|EPI_ISL_2333101\|EPI_ISL_2333102\|EPI_ISL_2333103\|EPI_ISL_2333104\|EPI_ISL_2333105\|EPI_ISL_2333106\|EPI_ISL_2333107\|EPI_ISL_2333108\|EPI_ISL_2333109\|EPI_ISL_2333110\|EPI_ISL_2333111\|EPI_ISL_2333112\|EPI_ISL_2333113\|EPI_ISL_2333114\|EPI_ISL_2333115\|EPI_ISL_2333116\|EPI_ISL_2333117\|EPI_ISL_2333178\|EPI_ISL_2333182\|EPI_ISL_2333183\|EPI_ISL_2333184\|EPI_ISL_2333185\|EPI_ISL_2333186\|EPI_ISL_2333187\|EPI_ISL_2333190\|EPI_ISL_2333191\|EPI_ISL_2333192\|EPI_ISL_2333193\|EPI_ISL_2333194\|EPI_ISL_2333195\|EPI_ISL_2333196\|EPI_ISL_2333202\|EPI_ISL_2333203\|EPI_ISL_2333214\|EPI_ISL_2333215\|EPI_ISL_2333216\|EPI_ISL_2333217\|EPI_ISL_2333516\|EPI_ISL_2333517\|EPI_ISL_2333518\|EPI_ISL_2333519\|EPI_ISL_2333520\|EPI_ISL_2333521\|EPI_ISL_2333522\|EPI_ISL_2333523\|EPI_ISL_2333524\|EPI_ISL_2333525\|EPI_ISL_2333526\|EPI_ISL_2333527\|EPI_ISL_2333528\|EPI_ISL_2357883" results/gisaid-$DATE_TREE-$X-metadata-sampled-1.tsv > results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv

	./scripts/phylogeo_sankoff_general --in results/gisaid-$DATE_TREE-$X-1.tree --out results/gisaid-$DATE_TREE-$X-2.tree --metadata results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv --location_label $STATE non$STATE --cond 2 "==" Germany  --print-annotation false --print-internal-node-label false --single-child false --samples results/gisaid-$DATE_TREE-$X-samples1.txt


	./scripts/contract_short_branch --in results/gisaid-$DATE_TREE-$X-2.tree --out results/gisaid-$DATE_TREE-$X-cont.tree --short 5e-6 --location_label $STATE non$STATE --print-annotation false --print-internal-node-label false --contract_leaf_enabled false

	./scripts/fill-template.py --template data/X.fixedRootPrior.skygrid-template-thorney.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples1.txt --tree results/gisaid-$DATE_TREE-$X-cont.tree --tree_init results/gisaid-$DATE_TREE-$X-2.tree --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.xml --loc Germany --date $DATE_TREE  --clock_rate ${SUB_TREE_CLOCK_RATE[$X]} --chain 300000000 # --cutoff 2.25 --grid_points $[ 225 * 52 / 100 ]
done

# COUNT number of samples from Germany/non-Germany
Rscript scripts/count-by-country.R Germany `for X in $SUB_TREES; do echo results/gisaid-$DATE_TREE-$X-samples1.txt results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv; done`
#for X in $SUB_TREES; do
#	./scripts/phylogeo_sankoff_general --in results/gisaid-$DATE_TREE-$X-1.tree --out temp/results/gisaid-$DATE_TREE-$X-2.tree --metadata results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv --location_label $STATE non$STATE --cond 2 "==" Germany  --print-annotation true --print-internal-node-label false --single-child false --samples temp/results/gisaid-$DATE_TREE-$X-samples1.txt
#done
#for X in $SUB_TREES; do
#        cat temp/results/gisaid-$DATE_TREE-$X-2.tree | sed 's/\([A-Z_0-9]*\[[^]]*\]\)/\n\1/g' | grep '^[A-Z_0-9]\+' | grep 'location=Germany'
#done | wc -l



for X in $SUB_TREES; do 
for i in {31..35}; do 
	qsub -cwd -N beast-$X-$i -M foroughmand@gmail.com -l h_vmem=20G,mem_free=20G,s_vmem=20G -pe smp 3 -o results/beast/run/out/$X-$i -j y run-beast-st.sh $DIR/results/beast/ run/$X-$i $DATE_TREE $X
	#JOB_ID=$RANDOM
	#sem -j 20 bash run-beast-st.sh $DIR/results/beast/ run/$X-$i $DATE_TREE $X 2>results/beast/run/out/$X-$i/beast-$X-$i.e$JOB_ID >results/beast/run/out/$X-$i/beast-$X-$i.o$JOB_ID
done
done

# Status
for X in $SUB_TREES; do 
for i in {31..35}; do 
GRID_ID=`ls results/beast/run/out/$X-$i/beast-$X-$i.e* | sed 's/.*e//' | sort | tail -n1`
LAST_LINE=`tail -n1 results/beast/run/out/$X-$i/beast-$X-$i.o$GRID_ID`
if [ "$LAST_LINE" != "done" ]; then
tail -n 5 results/beast/run/out/$X-$i/beast-$X-$i.e$GRID_ID results/beast/run/out/$X-$i/beast-$X-$i.o$GRID_ID
fi
done
done

#show status of running mcmc jobs
STAT_MCMC_RUN_DONE=""
STAT_MCMC_RUN_NO=""
#for i in {1..5}; do for X in B.1.1.70 B.1.1.317; do 
#for i in {1..10}; do for X in $SUB_TREES; do 
for i in {31..35}; do for X in $SUB_TREES; do 
if [[ -f results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-20210602.trees ]]; then
tail -n 1 results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-20210602.trees | grep -q '^End;'
if [ $? == 1 ]; then
echo -n "$X-$i ";
tail -n 1 results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-20210602.log | cut -f 1
else
STAT_MCMC_RUN_DONE="$STAT_MCMC_RUN_DONE $X-$i ";
echo -n "D: $X-$i "
tail -n 1 results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-20210602.log | cut -f 1
fi 
else
STAT_MCMC_RUN_NO="$STAT_MCMC_RUN_NO $X-$i";
fi
done done
echo "Done: $STAT_MCMC_RUN_DONE"
echo "Not started: $STAT_MCMC_RUN_NO"


for X in $SUB_TREES; do 
#for i in {1..5}; do 
for i in {31..35}; do 
loganalyser results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.log
done
done

#combining and making a log file for each subtree
#X=B.1.1.7
for X in $SUB_TREES; do 
logcombiner `
for i in {31..35}; do
cp  results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.log results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-fixed.log
#vim results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-fixed.log
#loganalyser results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-fixed.log
echo results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-fixed.log
done` results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.log >/dev/null
done


for X in $SUB_TREES; do 
MCMC_FRP_LOG_FILES=""
RANGE=`seq 31 35`
for i in $RANGE; do 
python scripts/resample.py --burnin 15000000 --rate 100000 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sampled.trees >/dev/null
python scripts/resample.py --burnin 15000000 --count 500 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sub500.trees >/dev/null
MCMC_FRP_LOG_FILES="$MCMC_FRP_LOG_FILES results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sub500.trees"
done
#mkdir results/beast/run/$X-all/
logcombiner -trees $MCMC_FRP_LOG_FILES results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.trees >/dev/null
done

#second round of MCMC, DTA

for X in $SUB_TREES; do 
	./scripts/fill-template.py --template data/X-DTA-template.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples1.txt --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X-DTA-$DATE_TREE.xml --loc Germany --date $DATE_TREE 
done

# run at grid.bifo
for X in $SUB_TREES; do 
for i in {31..32}; do 
qsub -cwd -N DTA-$X-$i -l h_vmem=64G,mem_free=20G,s_vmem=20G -pe smp 5 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA-st.sh run/$X-$i $DATE_TREE $X
#PID=$RANDOM$RANDOM
#sem -j 20 bash run-beast-DTA-st.sh run/$X-$i $DATE_TREE $X 2>results/beast/run/out/$X-$i/DTA-$X-$i.e$PID >results/beast/run/out/$X-$i/DTA-$X-$i.o$PID 
done
done


mkdir results/beast/run/all/
for X in $SUB_TREES; do 
MCMC_DTA_LOG_FILES=""
MCMC_DTA_TREE_FILES=""
#for i in {1..2}; do 
for i in {31..32}; do 
python scripts/resample.py --burnin 500000 --rate 4500 --tree results/beast/run/$X-$i/$X-DTA-$DATE_TREE.trees --out results/beast/run/$X-$i/$X-DTA-$DATE_TREE.sub4500.trees
MCMC_DTA_LOG_FILES="$MCMC_DTA_LOG_FILES results/beast/run/$X-$i/$X-DTA-$DATE_TREE.log"
MCMC_DTA_TREE_FILES="$MCMC_DTA_TREE_FILES results/beast/run/$X-$i/$X-DTA-$DATE_TREE.sub4500.trees"
done
echo logcombiner -trees $MCMC_DTA_TREE_FILES results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees
logcombiner -trees $MCMC_DTA_TREE_FILES results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees
echo logcombiner $MCMC_DTA_LOG_FILES results/beast/run/all/$X-DTA-$DATE_TREE.combined.log
logcombiner $MCMC_DTA_LOG_FILES results/beast/run/all/$X-DTA-$DATE_TREE.combined.log
done


mkdir results/beast/run/all
for i in {1..1}; do 
for X in $SUB_TREES; do 
#treeannotator -type mcc results/beast/run/$X-$i/$X-DTA-$DATE_TREE.sub4500.trees results/beast/run/all/$X-DTA-$DATE_TREE.MCC.trees
#$i is useless, should be removed
qsub -cwd -N ta-$X-$i -l h_vmem=64G,mem_free=20G,s_vmem=20G -pe smp 2 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-treeannotator-st.sh $DATE_TREE $X $i
#bash run-treeannotator-st.sh $DATE_TREE $X $i 
done
done

for X in $SUB_TREES; do 
xz -f -k -v -T0 results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees
xz -f -v -T0 results/beast/run/all/$X-DTA-$DATE_TREE.combined.log
done

R
cd reports
rmarkdown::render('extractLineages.Rmd');
rmarkdown::render('extractLineagesMCC.Rmd');


####
# BEGIN: The onetree model: (based on epa-ng injection of unsampled samples)
#   this is not currently the main method.
####

#unsampled: create fasta files
#SUB_TREES="B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.351 A C B.1.617.2"
scripts/extract-unsampled --location_label $STATE non$STATE --in results/gisaid-$DATE_TREE-$STATE.tree --samples results/gisaid-$DATE_TREE-samples.txt --unsampled results/gisaid-$DATE_TREE-unsampled.txt $SUB_TREES --cond 2 "==" Germany --metadata results/gisaid-$DATE_METADATA-metadata.tsv
cat results/gisaid-$DATE_TREE-unsampled.txt | grep -v "NA$" > results/gisaid-$DATE_TREE-unsampled-subtree.txt
./scripts/extract-metadata.R results/gisaid-$DATE_TREE-unsampled-subtree.txt results/gisaid-$DATE_METADATA-metadata.tsv results/gisaid-$DATE_TREE-metadata-unsampled.tsv
cat results/gisaid-$DATE_TREE-metadata-sampled.tsv results/gisaid-$DATE_TREE-metadata-unsampled.tsv  > results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv

##### TESTIGN ......
#####   ADDIG UNSAMPLED to the MCC tree
# tree  : results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree
# seq   : results/gisaid-$DATE_TREE-$X-seq0.fa
# seq.aln: results/fasttree/$DATE_TREE-$X/gisaid-$DATE_TREE-ft.xz
# unseq : results/gisaid-$DATE_TREE-$X-unseq0.fa
# 

#for X in $SUB_TREES; do
#echo MAFFT $X
#mafft --thread 26 --6merpair --addfragments results/gisaid-$DATE_TREE-$X-unseq0.fa --keeplength <(xzcat results/fasttree/$DATE_TREE-$X/gisaid-$DATE_TREE-ft.xz) > results/beast/unsampled/gisaid-$DATE_TREE-$X-unseq-seq.aln
#done

cut -f3 results/gisaid-$DATE_TREE-metadata-unsampled.tsv | tail -n +2  > results/mmsa-$DATE_TREE-unsample-ids.txt 
for X in $SUB_TREES; do cat results/gisaid-$DATE_TREE-$X-samples0.txt; done > results/gisaid-$DATE_TREE-all-samples0.txt
./scripts/alignment-filter ../../../data/data/mmsa_20210622_masked.fa <(cat results/mmsa-$DATE_TREE-unsample-ids.txt results/gisaid-$DATE_TREE-all-samples0.txt) > data/mmsa-$DATE_SEQ-sampled-unsampled0.fa 
mkdir results/beast/unsampled/

mkdir results/beast/run/lin-ius/ results/beast/run/lin-ius/out/
mkdir results/beast/run/lin-ius/aln/ results/beast/run/lin-ius/alnq/
for X in $SUB_TREES; do
	mkdir -p results/beast/unsampled/alignment-$DATE_TREE-$X-epa/
	echo "Enriching tree $X"


	./scripts/alignment-filter data/mmsa-$DATE_SEQ-sampled-unsampled0.fa <( grep "$X$" results/gisaid-$DATE_TREE-unsampled-subtree.txt | cut -f 1 ) > results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa 
	./scripts/alignment-filter data/mmsa-$DATE_SEQ-sampled-unsampled0.fa results/gisaid-$DATE_TREE-$X-samples0.txt > results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa 

	csplit -q --prefix=results/beast/run/lin-ius/aln/alignment-$DATE_TREE-$X-epa-reference- --suffix-format=%06d.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa '/^>/' '{*}'
	rm results/beast/run/lin-ius/aln/alignment-$DATE_TREE-$X-epa-reference-000000.fasta
	csplit -q --prefix=results/beast/run/lin-ius/alnq/alignment-$DATE_TREE-$X-epa-query- --suffix-format=%06d.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa '/^>/' '{*}'
	rm results/beast/run/lin-ius/alnq/alignment-$DATE_TREE-$X-epa-query-000000.fasta
#ls -S results/beast/run/lin-ius/alnq/ | grep "alignment-$DATE_TREE-$X-epa-query-" | 
ls -S results/beast/run/lin-ius/aln/ | grep "alignment-$DATE_TREE-$X-epa-reference-" | awk '{print "'results/beast/run/lin-ius/aln/'" $0}' | grep -v '\-000000.fasta' > results/beast/run/lin-ius/out/alignment-files-$X-reference.txt
	ls -S results/beast/run/lin-ius/alnq/ | grep "alignment-$DATE_TREE-$X-epa-query-" | awk '{print "'results/beast/run/lin-ius/alnq/'" $0}' | grep -v '\-000000.fasta' > results/beast/run/lin-ius/out/alignment-files-$X-query.txt
	mash sketch -o results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-reference.fasta results/beast/run/lin-ius/aln/alignment-$DATE_TREE-$X-epa-reference-*.fasta  2>>e
	mash sketch -o results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta -l results/beast/run/lin-ius/out/alignment-files-$X-query.txt 2>>e
	#mash dist -v 0.05 -p 52 results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-reference.fasta.msh results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt	
	mash dist -v 0.05 -p 120 results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-reference.fasta.msh results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt.gz
	#zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt.gz | python scripts/convert-filename-to-id.py | gzip > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz 
	zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt.gz | scripts/convert-filename-to-id | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz 
	#calculate distances between sequences
	#Rscript scripts/filter-by-country.R "$STATE" results/gisaid-$DATE_TREE-metadata-unsampled.tsv results/beast/unsampled/unsampled-$X-in.txt
	#Calculating distances between unsampled Germany sequences
	#Following should be executed for ius-3-aux
	#mash dist -v 0.05 -p 120 results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X.txt.gz
	#zcat results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X.txt.gz | scripts/convert-filename-to-id | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X-id.txt.gz
done

for X in $SUB_TREES; do
	echo "tree mash 2: $X"
	#Rscript scripts/filter-by-country.R "$STATE" results/gisaid-$DATE_TREE-metadata-unsampled.tsv results/beast/unsampled/unsampled-$X-in.txt
	#Calculating distances between unsampled Germany sequences
	#Following should be executed for ius-3-aux
	mash dist -v 0.05 -p 120 results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X.txt.gz
	zcat results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X.txt.gz | scripts/convert-filename-to-id | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X-id.txt.gz

	#qsub -cwd -N mash2-$X -l h_vmem=10G,mem_free=10G,s_vmem=10G -pe smp 3 -o results/beast/run/out/$X-1 -e results/beast/run/out/$X-1 run-mash2.sh $X

done

## fast extractLineage from sequences (add lineages to normal trees):
# INPUT:
#   tree: results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree
#   tree sample sequences: results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa
#   query seqeunces: results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa
#   results/gisaid-$DATE_TREE-metadata-sampled.tsv 
#   results/gisaid-$DATE_TREE-metadata-unsampled.tsv 
#   results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv
#   output folder: 
#     clusterSamples_DTA_MCC_NA.tsv
#     clusters_DTA_MCC_NA.tsv
# OUTPUT: results/beast/run/lin-ius/ (injected unsampled)
# we assume that input files are already created via incomplete run of epa-ng (previous method for lin-rich)
SUB_TREES="A B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.258 B.1.351 C"

#why some file names are different?
for X in $SUB_TREES; do
	mv results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fasta
	mv results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fasta
	mv results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa
	mv results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa
done
#for X in B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.351 A C B.1.617.2; do
for X in $SUB_TREES; do
	echo "Tree $X "
	cp results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus
	./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation false --single-child false --ilabel true --set-tip-location false
	./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE-rich.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation true --single-child false --ilabel true --set-tip-location false
	
	sed 's/Germany+nonGermany/nonGermany/g' results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus > results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus
	scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles_0.5.tsv

	# Nothing except sampled sequences
	#   Output folder: results/beast/run/lin-ius-no/
	#scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( cat /dev/null ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius-no/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius-no/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius-no/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius-no/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius-no/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius-no/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius-no/out/$X-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage true

	# For each sample, most similar among Germany and nonGermany is calculated, if it is from nonGermany, they make a new lineage
	#   Output folder: results/beast/run/lin-ius-2/
	#mkdir -p results/beast/run/lin-ius-2/out/
	#scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius-2/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius-2/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius-2/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius-2/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius-2/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius-2/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius-2/out/$X-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage true

	# Only allow identicals to be matched
	#   Output folder: results/beast/run/lin-ius-3/
	mkdir -p results/beast/run/lin-ius-3/out/
	scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius-3/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius-3/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage false --dist-threshold 0.0
	#between samples mkdir -p results/beast/run/lin-ius-3-aux/out/
	#scripts/lineage-importation-inject-unsampled --no-tree true --dist <( zcat results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius-3/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius-3-aux/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius-3-aux/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius-3-aux/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius-3-aux/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3-aux/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3-aux/out/$X-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage false --dist-threshold 0.0 --filter-out <( cut -f3 results/beast/run/lin-ius-3/out/$X-clusterSamples_DTA_MCC_NA.tsv ) 
	#


	#calculate lineage entrance to Germany, 
	#attach each new node to the lineage with most similar sample with tmrca > sample date
	#file is generated with sample-id -> cluster name
	#R create outputs

	#Calculating lineages without injection of unsampled
	#scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( echo -n "" ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-mcc/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-mcc/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-mcc/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-mcc/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-mcc/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-mcc/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-mcc/out/$X-clusters_DTA_MCC_singles_0.5.tsv
done

##Statistics for manuscript
#mkdir -p temp/lin-ius-3/out/
#for X in $SUB_TREES; do
#	scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder temp/lin-ius-3/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters temp/lin-ius-3/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples temp/lin-ius-3/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single temp/lin-ius-3/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 temp/lin-ius-3/out/$X-clusters_DTA_MCC_0.5.tsv temp/lin-ius-3/out/$X-clusterSamples_DTA_MCC_0.5.tsv temp/lin-ius-3/out/$X-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage false --dist-threshold 0.0 --log temp/lin-ius-3/$X.log &
#done

function lin_ius_clusters_combine {
LINIUS_FOLDER=$1
(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/$LINIUS_FOLDER/out/$X-clusters_DTA_MCC_NA.tsv
for X in $SUB_TREES; do
	tail -n +2 results/beast/run/$LINIUS_FOLDER/out/$X-clusters_DTA_MCC_NA.tsv 
done ) > results/beast/run/$LINIUS_FOLDER/clusters_DTA_MCC_NA.tsv

(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/$LINIUS_FOLDER/out/$X-clusterSamples_DTA_MCC_NA.tsv
for X in $SUB_TREES; do
	tail -n +2 results/beast/run/$LINIUS_FOLDER/out/$X-clusterSamples_DTA_MCC_NA.tsv
done ) > results/beast/run/$LINIUS_FOLDER/clusterSamples_DTA_MCC_NA.tsv

(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/$LINIUS_FOLDER/out/$X-clusters_DTA_MCC_0.5.tsv
for X in $SUB_TREES; do
	tail -n +2 results/beast/run/$LINIUS_FOLDER/out/$X-clusters_DTA_MCC_0.5.tsv 
done ) > results/beast/run/$LINIUS_FOLDER/clusters_DTA_MCC_0.5.tsv

(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/$LINIUS_FOLDER/out/$X-clusterSamples_DTA_MCC_0.5.tsv
for X in $SUB_TREES; do
	tail -n +2 results/beast/run/$LINIUS_FOLDER/out/$X-clusterSamples_DTA_MCC_0.5.tsv
done ) > results/beast/run/$LINIUS_FOLDER/clusterSamples_DTA_MCC_0.5.tsv

}

lin_ius_clusters_combine lin-ius
#lin_ius_clusters_combine lin-ius-2
lin_ius_clusters_combine lin-ius-3
#lin_ius_clusters_combine lin-ius-no
#lin_ius_clusters_combine lin-ius-3-aux

git add results/beast/run/lin-ius/cluster*
git add results/beast/run/lin-ius-3/cluster*
git add data/*.csv 
git add results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv

# For Denis' Team
#STATISTICS:
LARGEST_NAMES=16:116:164:181:237

mkdir -p temp/lin-ius-3/out
for X in $SUB_TREES; do
scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder temp/lin-ius-3/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters temp/lin-ius-3/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples temp/lin-ius-3/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single temp/lin-ius-3/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 temp/lin-ius-3/out/$X-clusters_DTA_MCC_0.5.tsv temp/lin-ius-3/out/$X-clusterSamples_DTA_MCC_0.5.tsv temp/lin-ius-3/out/$X-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage false --dist-threshold 0.0 --log temp/lin-ius-3/$X-2.log --do-assign true
done

grep 'lin assigned' `for X in $SUB_TREES; do echo "temp/lin-ius-3/$X-2.log"; done` | awk '{print $4}'


X=B.1.1.7
LARGEST_REXP=`echo $LARGEST_NAMES | sed 's/\([0-9]\+\)/DTA_MCC_\1=/g' | sed 's/:/\\\\|/g'`
grep "$LARGEST_REXP" temp/lin-ius-3/$X-2.log | sed 's/$/:1.00;/' > temp/lin-ius-3/largest.txt
split -l 5 -a 1 temp/lin-ius-3/largest.txt temp/lin-ius-3/largest-
mv temp/lin-ius-3/largest-a temp/lin-ius-3/largest-a.txt
mv temp/lin-ius-3/largest-b temp/lin-ius-3/largest-b.txt 

cat results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa > temp/lin-ius-3/$X-seqs.fa
cat temp/lin-ius-3/largest-ids.txt | while read a; do
N=`echo $a | cut -d'=' -f1`
echo $a | cut -d' ' -f2- | sed 's/ /\n/g' | grep -v '^$' > temp/lin-ius-3/$N-ids.txt
./scripts/alignment-filter temp/lin-ius-3/$X-seqs.fa temp/lin-ius-3/$N-ids.txt > temp/lin-ius-3/$N-seqs.fa
done

#BUILD DISTANCE MATRIX BETWEEN SAMPLED AND ADDED UNSAMPLED SEQUENCES
#for X in $SUB_TREES; do
#cat temp/lin-ius-3/$X-2.log | sed -n '/query_min_dist_lineage_half/,$p' | sed '/lin assigned/q' | sed '1d' | sed '$d' 
#done | cut -d' ' -f3 | sort | uniq > temp/lin-ius-3/query-ids.txt
##
#for X in $SUB_TREES; do
#cat temp/lin-ius-3/$X-2.log | sed -n '/query_min_dist_lineage_half/,$p' | sed '/lin assigned/q' | sed '1d' | sed '$d' 
#done | cut -d' ' -f4 | sort | uniq > temp/lin-ius-3/sample-idx.txt
#
for X in $SUB_TREES; do
cat temp/lin-ius-3/$X-2.log | sed -n '/query_min_dist_lineage_half/,$p' | sed '/lin assigned/q' | sed '1d' | sed '$d'
done > temp/lin-ius-3/match-mash-info.txt

for X in $SUB_TREES; do
cat temp/lin-ius-3/$X-2.log | sed -n '/query_min_dist_lineage_half/,$p' | sed '/lin assigned/q' | sed '1d' | sed '$d' 
done | cut -d' ' -f3,4  > temp/lin-ius-3/match-ids.txt

./scripts/alignment-filter data/mmsa-$DATE_SEQ-sampled-unsampled0.fa <( cat temp/lin-ius-3/query-ids.txt temp/lin-ius-3/sample-idx.txt ) > temp/lin-ius-3/query-sample-seqs.fa 

srun python scripts/compare-sequences.py temp/lin-ius-3/query-sample-seqs.fa temp/lin-ius-3/match-ids.txt > temp/lin-ius-3/sample-injected-distance.txt 2>temp/lin-ius-3/sample-injected-distance.err


for X in $SUB_TREES; do
	scripts/filter-alignment-distance <( cat temp/lin-ius-3/query-ids.txt temp/lin-ius-3/sample-idx.txt ) <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz )
done > temp/lin-ius-3/alignment-distance-injected.txt
python scripts/compare-sequences-mash.py temp/lin-ius-3/alignment-distance-injected.txt temp/lin-ius-3/match-mash-info.txt > temp/lin-ius-3/sample-injected-distance-mash.txt 2> temp/lin-ius-3/sample-injected-distance-mash.err



#CALCULATE number of subsampled sequences inside and outside Germany
for X in $SUB_TREES; do
	grep '^tree' results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus | sed 's/\([0-9]*\[[^]]*\]\)/\n\1/g' | grep '^[0-9]\+' | grep 'location="nonGermany"'
done | wc -l
for X in $SUB_TREES; do
	grep '^tree' results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus | sed 's/\([0-9]*\[[^]]*\]\)/\n\1/g' | grep '^[0-9]\+' | grep 'location="Germany"'
done | wc -l


#######################################
## NEW DTA for more than two states  ##
#######################################

mkdir results/dtamulti/

for X in $SUB_TREES; do 
./scripts/metadata-remove-na.R results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv results/dtamulti/gisaid-$DATE_TREE-$X-metadata-sampled-no-na.tsv
./scripts/phylogeo_sankoff_general --in  results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --save-nexus --out results/dtamulti/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --metadata results/dtamulti/gisaid-$DATE_TREE-$X-metadata-sampled-no-na.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany  --print-annotation false --print-internal-node-label false --single-child false --nosankoff --remove-invalid-children true
done

./scripts/merge-trees-for-multidta --trees `for X in $SUB_TREES; do echo results/dtamulti/$X-DTA-$DATE_TREE.MCC.tree.2.nexus; done` --root-date 1900 --out results/dtamulti/sampled.fixedRootPrior.skygrid-$DATE_TREE.trees -l \"Germany\" \"nonGermany\" --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --samples results/dtamulti/sampled.fixedRootPrior.skygrid-$DATE_TREE.samples

./scripts/fill-template.py --template data/X-DTAmulti-template.xml --name sampled --sample results/dtamulti/sampled.fixedRootPrior.skygrid-$DATE_TREE.samples --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/dtamulti/sampled-DTAmulti-$DATE_TREE.xml --loc Germany --date $DATE_TREE --loc_depth

# run at grid.bifo
for i in {31..32}; do 
sbatch --job-name=mdta-$i run-beast-DTAmulti.sh run/sampled-$i $DATE_TREE sampled
done

mkdir results/dtamulti/run/all results/dtamulti/run/lin-ius-3/
MCMC_DTA_LOG_FILES=""
MCMC_DTA_TREE_FILES=""
for i in {31..32}; do 
srun python scripts/resample.py --burnin 500000 --rate 4500 --tree results/dtamulti/run/sampled-$i/sampled-DTA-$DATE_TREE.trees --out results/dtamulti/run/sampled-$i/sampled-DTA-$DATE_TREE.sub4500.trees
MCMC_DTA_LOG_FILES="$MCMC_DTA_LOG_FILES results/dtamulti/run/sampled-$i/sampled-DTA-$DATE_TREE.log"
MCMC_DTA_TREE_FILES="$MCMC_DTA_TREE_FILES results/dtamulti/run/sampled-$i/sampled-DTA-$DATE_TREE.sub4500.trees"
done
echo logcombiner -trees $MCMC_DTA_TREE_FILES results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.trees
srun logcombiner -trees $MCMC_DTA_TREE_FILES results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.trees
echo logcombiner $MCMC_DTA_LOG_FILES results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.log
srun logcombiner $MCMC_DTA_LOG_FILES results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.log

loganalyser results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.log > results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.log.analyzed


alias srun='srun -w ariel --cpus-per-task=2 --qos verylong --time=24:00:00 '

srun --mem=100G --cpus-per-task=2 --qos verylong --time=24:00:00 ../../bin/beast/bin/treeannotator -type mcc results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.trees  results/dtamulti/run/all/sampled-DTA-$DATE_TREE.MCC.nexus

srun ./scripts/lineage-importation-extract-multidta --tree results/dtamulti/run/all/sampled-DTA-$DATE_TREE.MCC.nexus --lineage-samples results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv --out results/dtamulti/run/lin-ius-3/clusterMovement_DTA_MCC_0.5.0.tsv -l \"Baden-Wurttemberg\" \"Bavaria\" \"Berlin\" \"Brandenburg\" \"Bremen\" \"Hamburg\" \"Hesse\" \"Lower Saxony\" \"Mecklenburg-Western Pomerania\" \"North Rhine-Westphalia\" \"Rhineland-Palatinate\" \"Saarland\" \"Saxony\" \"Saxony-Anhalt\" \"Schleswig-Holstein\" \"Thuringia\" --root-date 1900
sed 's/"//g' < results/dtamulti/run/lin-ius-3/clusterMovement_DTA_MCC_0.5.0.tsv > results/dtamulti/run/lin-ius-3/clusterMovement_DTA_MCC_0.5.tsv

##  POSTERIOR DISTRIBUTION SAMPLES OF MULTIDTA
rm -rf results/dtamulti/run/lin-ius-3e
mkdir -p results/dtamulti/run/lin-ius-3e/sample-tree/ results/dtamulti/run/lin-ius-3e/movement/ results/dtamulti/run/lin-ius-3e/log/

alias srun='srun --cpus-per-task=2 --qos verylong --time=24:00:00 '

srun python scripts/tree-separate.py results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.trees results/dtamulti/run/lin-ius-3e/sample-tree/tree-NAME.single.nexus

for TREE_FILE in results/dtamulti/run/lin-ius-3e/sample-tree/*.single.nexus ; do
TREE_NAME=`echo $(basename $TREE_FILE) | sed 's/\.single\.nexus//'`
echo $TREE_NAME $TREE_FILE
#sbatch --mem=50g --time=1:00:00 --cpus-per-task 1 -e results/dtamulti/run/lin-ius-3e/log/eff-$TREE_NAME.log -o results/dtamulti/run/lin-ius-3e/log/eff-$TREE_NAME.log scripts/generate-effectiveness.sh $DATE_TREE $STATE $TREE_NAME
sbatch --mem=50g --time=1:00:00 --cpus-per-task 1 -e results/dtamulti/run/lin-ius-3e/log/eff-$TREE_NAME.err.log -o results/dtamulti/run/lin-ius-3e/log/eff-$TREE_NAME.out.log run-2-movement.sh $TREE_NAME
done

TREE_NAME=tree-997416000
head -n 1 results/dtamulti/run/lin-ius-3e/movement/$TREE_NAME.clusterMovement_DTA_MCC_0.5.0.tsv | sed 's/$/\ttree/' > results/dtamulti/run/lin-ius-3e/all.clusterMovement_DTA_MCC_0.5.0.tsv

for TREE_FILE in results/dtamulti/run/lin-ius-3e/sample-tree/*.single.nexus ; do
TREE_NAME=`echo $(basename $TREE_FILE) | sed 's/\.single\.nexus//'`
tail -n+2 results/dtamulti/run/lin-ius-3e/movement/$TREE_NAME.clusterMovement_DTA_MCC_0.5.0.tsv | sed 's/$/\t'$TREE_NAME'/' >> results/dtamulti/run/lin-ius-3e/all.clusterMovement_DTA_MCC_0.5.0.tsv
done

## POSTERIOR DISTRIBUTION SAMPLES OF MULTIDTA - STARTING FROM SAMPLES (NOT MCC) OF DTA

##  POSTERIOR DISTRIBUTION SAMPLES OF DTA

rm -rf results/beast/run/lin-ius-3e/

#extract trees results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees to results/beast/run/lin-ius-3e/sample-tree
mkdir -p results/beast/run/lin-ius-3e/out results/beast/run/lin-ius-3e/sample-tree results/beast/run/lin-ius-3e/log results/beast/run/lin-ius-3e/mcc-tree/  results/beast/run/lin-ius-3e/out-tree/ 
for X in $SUB_TREES; do
	echo $X ...
	srun python scripts/tree-separate.py results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees results/beast/run/lin-ius-3e/sample-tree/$X-NAME.nexus
	#cat results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees | sed -n '/^tree /q;p' > results/beast/run/lin-ius-3e/sample-tree/$X-pre-trees.tmp
	#cat results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees | while read -r a; do
	#	if [[ "$a" =~ ^tree\  ]]; then
	#		name=$X-$(echo "$a" | awk '{print $2}' | sed 's/STATE_//')
	#		echo $name
	#		(
	#		cat results/beast/run/lin-ius-3e/sample-tree/$X-pre-trees.tmp
	#		echo $a
	#		echo "End;"
	#		) > results/beast/run/lin-ius-3e/sample-tree/$name.nexus
	#	fi
	#done
	#rm results/beast/run/lin-ius-3e/sample-tree/$X-pre-trees.tmp
done

for TREE_FILE in results/beast/run/lin-ius-3e/sample-tree/*.nexus ; do
TREE_NAME=`echo $(basename $TREE_FILE) | sed 's/\.nexus//'`
echo $TREE_NAME $TREE_FILE
sbatch --mem=50g --time=1:00:00 --cpus-per-task 1 -e results/beast/run/lin-ius-3e/log/eff-$TREE_NAME.log -o results/beast/run/lin-ius-3e/log/eff-$TREE_NAME.log scripts/generate-effectiveness.sh $DATE_TREE $STATE $TREE_NAME
done
# results in results/beast/run/lin-ius-3e/out/$TREE_NAME-effectiveness.tsv
library(dplyr)
library(tidyr)
library(bfp)
d <- read.table('results/beast/run/lin-ius-3e/all.tsv', sep='\t', header=TRUE, stringsAsFactor=FALSE)
colnames(d) <- c('name', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8',  'N9', 'N10', 'N11', 'N12', 'tree')
d<- d %>% filter(name != 'name')
ds <- d %>% mutate(name = ifelse(name == " no NRW,BV", "no NRW,BV", name)) %>% pivot_longer(num_range("N", 1:12), names_to='npi', values_to='eff') %>% mutate(eff = as.numeric(eff)) %>% group_by(name, npi) %>% summarise(eff.avg = mean(eff), eff.sd=sd(eff), hpd.l = empiricalHpd(eff, level=0.95)[1], hpd.u=empiricalHpd(eff, level=0.95)[2]) %>% pivot_wider(names_from=npi, values_from=c(eff.avg, eff.sd, hpd.l, hpd.u)) %>% select(sort(tidyselect::peek_vars())) %>% relocate(name) 
write.table(ds, 'results/beast/run/lin-ius-3e/summ.tsv', sep='\t', quote=FALSE, row.names=FALSE)


TREE_NAME=A-1000440000
head -n 1 results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_0.5.tsv | sed 's/$/\ttree/' > results/beast/run/lin-ius-3e/out/all-clusters_DTA_MCC_0.5.tsv 
head -n 1 results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_0.5.tsv | sed 's/$/\ttree/' > results/beast/run/lin-ius-3e/all-clusters_DTA_MCC_0.5.tsv
head -n 1 results/beast/run/lin-ius-3e/out/$TREE_NAME-clusterSamples_DTA_MCC_0.5.tsv| sed 's/$/\ttree/' > results/beast/run/lin-ius-3e/out/all-clusterSamples_DTA_MCC_0.5.tsv
head -n 1 results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_singles_0.5.tsv | sed 's/$/\ttree/' >> results/beast/run/lin-ius-3e/out/all-clusters_DTA_MCC_singles_0.5.tsv

for TREE_FILE in results/beast/run/lin-ius-3e/sample-tree/*.nexus ; do
TREE_NAME=`echo $(basename $TREE_FILE) | sed 's/\.nexus//'`
#echo $TREE_NAME $TREE_FILE
tail -n+2 results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_0.5.tsv  | sed 's/$/\t'$TREE_NAME'/' >> results/beast/run/lin-ius-3e/out/all-clusters_DTA_MCC_0.5.tsv 
tail -n+2 results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_0.5.tsv | sed 's/$/\t'$TREE_NAME'/' >> results/beast/run/lin-ius-3e/all-clusters_DTA_MCC_0.5.tsv
tail -n+2 results/beast/run/lin-ius-3e/out/$TREE_NAME-clusterSamples_DTA_MCC_0.5.tsv | sed 's/$/\t'$TREE_NAME'/' >> results/beast/run/lin-ius-3e/out/all-clusterSamples_DTA_MCC_0.5.tsv
tail -n+2 results/beast/run/lin-ius-3e/out/$TREE_NAME-clusters_DTA_MCC_singles_0.5.tsv | sed 's/$/\t'$TREE_NAME'/' >> results/beast/run/lin-ius-3e/out/all-clusters_DTA_MCC_singles_0.5.tsv
done

d <- read.table('results/beast/run/lin-ius-3e/out/all-clusterSamples_DTA_MCC_0.5.tsv', sep='\t', header=TRUE) %>% select(cluster, tree)
ds <- d %>% group_by(cluster, tree) %>% summarise(n=n())
dss <- ds %>% group_by(tree) %>% summarise(n_cluster=n())
ds %>% mutate(posterior_sample_id=str_replace(tree, '.*-', '')) %>% group_by(posterior_sample_id) %>% summarise(cluster=n()) %>% arrange(-cluster) %>% hist()



#for X in $SUB_TREES; do
#scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( cat /dev/null ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius-3e/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/run/lin-ius-3e/out-tree/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius-3e/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius-3e/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius-3e/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius-3e/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3e/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3e/out/$X-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage false --dist-threshold 0.0 --log results/beast/run/lin-ius-3e/$X-2.log --do-assign true
#done

