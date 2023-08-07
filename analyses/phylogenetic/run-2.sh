#DATE_TREE=$1
#DATE_METADATA=$2
#STATE=Germany

export LD_LIBRARY_PATH="/home/hforoughmand/miniconda3/envs/covid-uk/lib/"
export PATH=/net/viral_genomics/covid-lineage/germany-lineage-dynamics/bin/beast/bin:$PATH
cd /net/viral_genomics/covid-lineage/huge-lineage-dynamics/analyses/phylogenetic
DATE_TREE=20210602
DATE_METADATA=20210602
DATE_SEQ=20210602
STATE=Germany
TWD=/net/sgi/viral_genomics/hadi/tmp/
SUB_TREES="A B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.258 B.1.351 C"
#SUB_TREES="B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.1.214 B.1.177 B.1.160 B.1.221 B.1.36 B.1.258 B.1.351 P A C"
#SUB_TREES="$SUB_TREES B.1.617.2"
#SUB_TREES="$SUB_TREES B.1.258" #ADDED Oct 21 2021

#V1: download tree from GISAID
#V2: build the tree
#
#cleaning
grep -A1000 -F 'hCoV-19/Wuhan/WH04/2020' ../../../data/data/gisaid-$DATE_SEQ-raw.fa > ../../../data/data/gisaid-$DATE_SEQ-Wu1.fa
#edit then ..., label: hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05
DIR=`pwd`
mkdir $TWD/gisaid-$DATE_SEQ/
cd $TWD/gisaid-$DATE_SEQ/
ln -s $DIR/../../../data/data/gisaid-$DATE_SEQ-raw.fa .
bash $DIR/sarscov2phylo/clean_gisaid.sh -i gisaid-$DATE_SEQ-raw.fa -o gisaid-$DATE_SEQ.fa -t 30 
#qsub qrun.sh
mv gisaid-$DATE_SEQ.fa $DIR/../../../data/data/

cd $DIR

## SUBSAMPLING: 
#V3 (CURRENT): build tree as pangolin tree: 
./scripts/tree-build-as-pangolin --out results/pangolin-$DATE_TREE-r0.tree --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv --seqs ../../../data/data/gisaid-$DATE_SEQ.fa  -l Germany nonGermany
# we use following for filtering non-good samaple (e.g. sample without complete date)
./scripts/phylogeo_sankoff_general --in results/pangolin-$DATE_TREE-r0.tree --out results/gisaid-$DATE_TREE-$STATE.tree --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv --location_label $STATE non$STATE --cond 2 "==" Germany 
#./scripts/contract_short_branch --in results/gisaid-$DATE_TREE-$STATE.tree --out results/gisaid-$DATE_TREE-cont.tree --short 5e-6 --location_label $STATE non$STATE --print-annotation true --print-internal-node-label true 

./scripts/sample-evenly --in results/gisaid-$DATE_TREE-$STATE.tree --out results/gisaid-$DATE_TREE-sampled.tree~ -l $STATE non$STATE --samples-out results/gisaid-$DATE_TREE-samples.txt~ --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv --bucket-size 100 25 --special 5 $SUB_TREES --keep 1000 10000 $SUB_TREES 
cp results/gisaid-$DATE_TREE-sampled.tree~ results/gisaid-$DATE_TREE-sampled.tree
cp results/gisaid-$DATE_TREE-samples.txt~ results/gisaid-$DATE_TREE-samples.txt

# We can check the statistics of the trees and subtrees with following commands. The SUB_TREES variable is calculated beased on following statistics.
./scripts/tree-stat --in results/gisaid-$DATE_TREE.tree -l $STATE non$STATE --large 1000 --depth 10
./scripts/tree-stat --in results/gisaid-$DATE_TREE-sampled.tree -l $STATE non$STATE --large 500 --depth 30

./scripts/extract-metadata.R results/gisaid-$DATE_TREE-samples.txt ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv results/gisaid-$DATE_TREE-metadata-sampled.tsv
cut -f1 results/gisaid-$DATE_TREE-metadata-sampled.tsv | tail -n +2  > results/gisaid-$DATE_TREE-sample-names.txt 
./scripts/alignment-filter ../../../data/data/gisaid-$DATE_SEQ.fa results/gisaid-$DATE_TREE-sample-names.txt > ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled0.fa 
grep '>' ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled0.fa | sort | uniq -c | grep -v '^ *1 ' | sort -nr 
./scripts/rename-alignment-ids results/gisaid-$DATE_TREE-metadata-sampled.tsv < ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled0.fa > ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled.fa
#./scripts/alignment-filter ../../data/phylogenetic/gisaid-$DATE_SEQ.fa.xz results/gisaid-$DATE_TREE-samples.txt > ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled.fa 

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
	./scripts/alignment-filter ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled.fa results/gisaid-$DATE_TREE-$X-samples0.txt > results/gisaid-$DATE_TREE-$X-seq0.fa
	cat ../../../data/data/gisaid-$DATE_SEQ-Wu1.fa >> results/gisaid-$DATE_TREE-$X-seq0.fa
	#cd results/fasttree/$DATE_TREE-$X/
	#bash $DIR/sarscov2phylo/global_tree_gisaid.sh -i gisaid-$DATE_TREE-$X-seq0.fa -o gisaid-$DATE_TREE-$X-tree -t 20
	#qsub -cwd -l h_vmem=10G,mem_free=10G,s_vmem=10G -pe smp 30 -o log -e log qrun.sh $DATE_SEQ $DATE_TREE $X
	#cd -
done


SUB_TREES=?????
#run on grid
for X in $SUB_TREES; do
	DIR=`pwd` ;
	rm -rf results/fasttree/$DATE_TREE-$X/ ;
	mkdir -p results/fasttree/$DATE_TREE-$X/ ;
	ln -s $DIR/results/gisaid-$DATE_TREE-$X-seq0.fa results/fasttree/$DATE_TREE-$X/gisaid-$DATE_TREE-$X-seq0.fa ;
	cp scripts/qrun-fasttree.sh results/fasttree/$DATE_TREE-$X/qrun.sh ;
	cd results/fasttree/$DATE_TREE-$X/ ;
	qsub -cwd -N sarscov2phylo-$X -M foroughmand@gmail.com -l h_vmem=30G,mem_free=30G,s_vmem=30G -pe smp 1 -o log -e log qrun.sh $DATE_SEQ $DATE_TREE $X ;
# qsub -cwd -l h_vmem=64G,mem_free=20G,s_vmem=20G -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA.sh run/$X-$i $DATE_TREE $X# qsub -cwd -l h_vmem=64G,mem_free=20G,s_vmem=20G -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA.sh run/$X-$i $DATE_TREE $X
# qsub -cwd -l h_vmem=64G,mem_free=20G,s_vmem=20G -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA.sh run/$X-$i $DATE_TREE $X
# -pe smp 10 
	cd - ;
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
# I changed initial values, but return to the old one. Now the clock.rate is fixed
#		[B.1.258]=0.000075
#		[B.1.1.70]=0.000075
#		[B.1.177]=0.000075

#remove case EPI_ISL_2333117 which has date 1905-06-26

for X in $SUB_TREES; do
	cp results/fasttree/$DATE_TREE-$X/ft_FBP.tree results/gisaid-$DATE_TREE-$X-1.tree ;
	./scripts/extract-metadata.R results/gisaid-$DATE_TREE-$X-samples0.txt results/gisaid-$DATE_TREE-metadata-sampled.tsv results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv

	# B.1.1.70
	grep -v EPI_ISL_2333117 results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv > x; cat x > results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv
	grep -v EPI_ISL_2333526 results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv > x; cat x > results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv

	#B.1.1.177
	grep -v EPI_ISL_2333528 results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv > x; cat x > results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv
	

	# B.1.258
	grep -v EPI_ISL_2333525 results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv > x; cat x > results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv

	# B.1.1.7
	grep -v EPI_ISL_2357883 results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv > x; cat x > results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv

	./scripts/phylogeo_sankoff_general --in results/gisaid-$DATE_TREE-$X-1.tree --out results/gisaid-$DATE_TREE-$X-2.tree --metadata results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv --location_label $STATE non$STATE --cond 2 "==" Germany  --print-annotation false --print-internal-node-label false --single-child false --samples results/gisaid-$DATE_TREE-$X-samples1.txt


	./scripts/contract_short_branch --in results/gisaid-$DATE_TREE-$X-2.tree --out results/gisaid-$DATE_TREE-$X-cont.tree --short 5e-6 --location_label $STATE non$STATE --print-annotation false --print-internal-node-label false --contract_leaf_enabled false

	./scripts/fill-template.py --template data/X.fixedRootPrior.skygrid-template-thorney.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples1.txt --tree results/gisaid-$DATE_TREE-$X-cont.tree --tree_init results/gisaid-$DATE_TREE-$X-2.tree --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.xml --loc Germany --date $DATE_TREE  --clock_rate ${SUB_TREE_CLOCK_RATE[$X]} --chain 300000000 # --cutoff 2.25 --grid_points $[ 225 * 52 / 100 ]
done

# run at grid.bifo
#SUB_TREES="B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.1.214 B.1.177 B.1.160 B.1.221 B.1.36 B.1.258 B.1.351 P A C"
SUB_TREES="B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.351 P A C"
#NOTE: B.1.1.7 -> i={6..10}
cd /net/viral_genomics/covid-lineage/huge-lineage-dynamics/analyses/phylogenetic/
for X in $SUB_TREES; do 
#for X in B.1.1.70 B.1.177; do
#current jobs number are 31..35, previously it was 1..5 for all and 6..20 for B.1.1.7
for i in {31..35}; do 
#if [ "$X" == 'B.1.1.7' ]; then
##qsub -cwd -N beast-$X-$i -M foroughmand@gmail.com -l h_vmem=10G,mem_free=10G,s_vmem=10G -pe smp 5 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast.sh run/$X-$i $DATE_TREE $X
#echo "E";
#else
qsub -cwd -N beast-$X-$i -M foroughmand@gmail.com -l h_vmem=10G,mem_free=10G,s_vmem=10G -pe smp 3 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast.sh run/$X-$i $DATE_TREE $X
#fi
done
done

## continue B.1.1.7 for more steps
#X=B.1.1.7
#for i in {11..15}; do 
#ii=$(( i-5 ))
#mkdir results/beast/run/$X-$i/
#cp results/beast/run/$X-$ii/$X.fixedRootPrior.skygrid-$DATE_TREE.state results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-last.state
## did I change the main xml?!
#./scripts/fill-template.py --template data/X.fixedRootPrior.skygrid-template-thorney.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples1.txt --tree results/gisaid-$DATE_TREE-$X-cont.tree --tree_init results/gisaid-$DATE_TREE-$X-2.tree --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.xml --loc Germany --date $DATE_TREE  --clock_rate ${SUB_TREE_CLOCK_RATE[$X]} --chain 600000000 
#done
#
#for i in {11..15}; do 
#qsub -cwd -N beast-$X-$i -M foroughmand@gmail.com -l h_vmem=10G,mem_free=10G,s_vmem=10G -pe smp 3 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-more.sh run/$X-$i $DATE_TREE $X
#done
#
#for i in {16..20}; do 
#qsub -cwd -N beast-$X-$i -M foroughmand@gmail.com -l h_vmem=10G,mem_free=10G,s_vmem=10G -pe smp 3 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast.sh run/$X-$i $DATE_TREE $X
#done


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
#for i in {11..20}; do
#cp  results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.log results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-fixed.log
#vim results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-fixed.log
#loganalyser results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-fixed.log
echo results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-fixed.log
done` results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.log
done


for X in $SUB_TREES; do 
MCMC_FRP_LOG_FILES=""
RANGE=`seq 31 35`
#RANGE=`seq 1 5`
#if [[ "$X" == "B.1.1.7" ]]; then 
##RANGE=`seq 6 10`
#RANGE=`seq 11 20`
#fi
for i in $RANGE; do 
python scripts/resample.py --burnin 15000000 --rate 100000 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sampled.trees
python scripts/resample.py --burnin 15000000 --count 500 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sub500.trees
MCMC_FRP_LOG_FILES="$MCMC_FRP_LOG_FILES results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sub500.trees"
done
#mkdir results/beast/run/$X-all/
logcombiner -trees $MCMC_FRP_LOG_FILES results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.trees
done

##resample for some remaining cases ...
##
#for X in B.1.1.70 B.1.1.317 B.1.177; do
#MCMC_FRP_LOG_FILES=""
#RANGE=`seq 1 5`
#for i in $RANGE; do
#echo $X $i
#python scripts/resample.py --burnin 15000000 --rate 100000 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sampled.trees
#python scripts/resample.py --burnin 15000000 --count 500 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sub500.trees
#MCMC_FRP_LOG_FILES="$MCMC_FRP_LOG_FILES results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sub500.trees"
#done
#logcombiner -trees $MCMC_FRP_LOG_FILES results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.trees
#
#resample
#
#done
#
#for X in B.1.1.70 B.1.1.317 B.1.177; do
#	./scripts/fill-template.py --template data/X-DTA-template.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples1.txt --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X-DTA-$DATE_TREE.xml --loc Germany --date $DATE_TREE 
#done
#
#for X in B.1.1.70 B.1.1.317 B.1.177; do
#for i in {31..32}; do 
#bash run-beast-DTA.sh run/$X-$i $DATE_TREE $X 2>results/beast/run/out/$X-$i/err >results/beast/run/out/$X-$i/out
#done
#done
#
#for X in B.1.1.70 B.1.1.317 B.1.177; do
#MCMC_DTA_LOG_FILES=""
#MCMC_DTA_TREE_FILES=""
##for i in {1..2}; do 
#for i in {31..32}; do 
#python scripts/resample.py --burnin 500000 --rate 4500 --tree results/beast/run/$X-$i/$X-DTA-$DATE_TREE.trees --out results/beast/run/$X-$i/$X-DTA-$DATE_TREE.sub4500.trees
#MCMC_DTA_LOG_FILES="$MCMC_DTA_LOG_FILES results/beast/run/$X-$i/$X-DTA-$DATE_TREE.log"
#MCMC_DTA_TREE_FILES="$MCMC_DTA_TREE_FILES results/beast/run/$X-$i/$X-DTA-$DATE_TREE.sub4500.trees"
#done
#echo logcombiner -trees $MCMC_DTA_TREE_FILES results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees
#logcombiner -trees $MCMC_DTA_TREE_FILES results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees
#echo logcombiner $MCMC_DTA_LOG_FILES results/beast/run/all/$X-DTA-$DATE_TREE.combined.log
#logcombiner $MCMC_DTA_LOG_FILES results/beast/run/all/$X-DTA-$DATE_TREE.combined.log
#done
#
#for i in {1..1}; do 
#for X in B.1.1.70 B.1.1.317 B.1.177; do
##treeannotator -type mcc results/beast/run/$X-$i/$X-DTA-$DATE_TREE.sub4500.trees results/beast/run/all/$X-DTA-$DATE_TREE.MCC.trees
##$i is useless, should be removed
#bash run-treeannotator.sh $DATE_TREE $X $i
#done
#done
#
#for X in B.1.1.70 B.1.1.317 B.1.177; do
#	#waiting for B.1.1.70 B.1.1.7 B.1.177 
#	#done C B.1.1.519 B.1.1.317 B.1.221 B.1.36 B.1.351 B.1.160 A
#	#for X in ; do
#	echo "Enriching tree $X"
#	cp results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus
#	./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation false --single-child false --ilabel true --set-tip-location false
#	./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE-rich.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation true --single-child false --ilabel true --set-tip-location false
#
#	echo "Tree $X "
#	#calculate distances between sequences
#	#Rscript scripts/filter-by-country.R "$STATE" results/gisaid-$DATE_TREE-metadata-unsampled.tsv results/beast/unsampled/unsampled-$X-in.txt
#	csplit -q --prefix=results/beast/run/lin-ius/aln/alignment-$DATE_TREE-$X-epa-reference- --suffix-format=%06d.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa '/^>/' '{*}'
#	rm results/beast/run/lin-ius/aln/alignment-$DATE_TREE-$X-epa-reference-000000.fasta
#	csplit -q --prefix=results/beast/run/lin-ius/alnq/alignment-$DATE_TREE-$X-epa-query- --suffix-format=%06d.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa '/^>/' '{*}'
#	rm results/beast/run/lin-ius/alnq/alignment-$DATE_TREE-$X-epa-query-000000.fasta
#	#ls -S results/beast/run/lin-ius/alnq/ | grep "alignment-$DATE_TREE-$X-epa-query-" | 
#	ls -S results/beast/run/lin-ius/aln/ | grep "alignment-$DATE_TREE-$X-epa-reference-" | awk '{print "'results/beast/run/lin-ius/aln/'" $0}' | grep -v '\-000000.fasta' > results/beast/run/lin-ius/out/alignment-files-$X-reference.txt
#	ls -S results/beast/run/lin-ius/alnq/ | grep "alignment-$DATE_TREE-$X-epa-query-" | awk '{print "'results/beast/run/lin-ius/alnq/'" $0}' | grep -v '\-000000.fasta' > results/beast/run/lin-ius/out/alignment-files-$X-query.txt
#	mash sketch -o results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-reference.fasta results/beast/run/lin-ius/aln/alignment-$DATE_TREE-$X-epa-reference-*.fasta  2>>e
#	mash sketch -o results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta -l results/beast/run/lin-ius/out/alignment-files-$X-query.txt 2>>e
#	#mash dist -v 0.05 -p 52 results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-reference.fasta.msh results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt	
#	time mash dist -v 0.05 -p 120 results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-reference.fasta.msh results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt.gz
#	#zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt.gz | python scripts/convert-filename-to-id.py | gzip > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz 
#	zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt.gz | scripts/convert-filename-to-id | gzip > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz 
#	
#	sed 's/Germany+nonGermany/nonGermany/g' results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus > results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus
#	scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles.tsv
#done
#
##END

#second round of MCMC, DTA

for X in $SUB_TREES; do 
	./scripts/fill-template.py --template data/X-DTA-template.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples1.txt --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X-DTA-$DATE_TREE.xml --loc Germany --date $DATE_TREE 
done

# run at grid.bifo
#DATE_TREE=20210428
cd /net/viral_genomics/covid-lineage/huge-lineage-dynamics/analyses/phylogenetic/
for X in $SUB_TREES; do 
#for i in {1..2}; do 
for i in {31..32}; do 
qsub -cwd -N DTA-$X-$i -l h_vmem=64G,mem_free=20G,s_vmem=20G -pe smp 5 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA.sh run/$X-$i $DATE_TREE $X
done
done


# Removed 
#  No sample from Germany: B.1.1.214
#  Only singletone: B.1.258, P
SUB_TREES="B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.351 A C"
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
qsub -cwd -N ta-$X-$i -l h_vmem=64G,mem_free=20G,s_vmem=20G -pe smp 2 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-treeannotator.sh $DATE_TREE $X $i
done
done

for X in $SUB_TREES; do 
xz -f -k -v -T0 results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees
xz -f -v -T0 results/beast/run/all/$X-DTA-$DATE_TREE.combined.log
done

#DONE: remove results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees



R
cd reports
rmarkdown::render('extractLineages.Rmd');
rmarkdown::render('extractLineagesMCC.Rmd');


####
# BEGIN: The onetree model: (based on epa-ng injection of unsampled samples)
#   this is not currently the main method.
####

#unsampled: create fasta files
SUB_TREES="B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.351 A C B.1.617.2"
scripts/extract-unsampled --location_label $STATE non$STATE --in results/gisaid-$DATE_TREE-$STATE.tree --samples results/gisaid-$DATE_TREE-samples.txt --unsampled results/gisaid-$DATE_TREE-unsampled.txt $SUB_TREES --cond 2 "==" Germany --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv
cat results/gisaid-$DATE_TREE-unsampled.txt | grep -v "NA$" > results/gisaid-$DATE_TREE-unsampled-subtree.txt
./scripts/extract-metadata.R results/gisaid-$DATE_TREE-unsampled-subtree.txt ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv results/gisaid-$DATE_TREE-metadata-unsampled.tsv
cat results/gisaid-$DATE_TREE-metadata-sampled.tsv results/gisaid-$DATE_TREE-metadata-unsampled.tsv  > results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv
#cut -f1 results/gisaid-$DATE_TREE-metadata-unsampled.tsv | tail -n +2  > results/gisaid-$DATE_TREE-unsample-names.txt 
#./scripts/alignment-filter ../../../data/data/gisaid-$DATE_SEQ.fa results/gisaid-$DATE_TREE-unsample-names.txt > ../../data/phylogenetic/gisaid-$DATE_SEQ-unsampled0.fa 
#grep '>' ../../data/phylogenetic/gisaid-$DATE_SEQ-unsampled0.fa | sort | uniq -c | grep -v '^ *1 ' | sort -nr 
#./scripts/rename-alignment-ids results/gisaid-$DATE_TREE-metadata-unsampled.tsv < ../../data/phylogenetic/gisaid-$DATE_SEQ-unsampled0.fa > ../../data/phylogenetic/gisaid-$DATE_SEQ-unsampled.fa
#for X in $SUB_TREES; do
#	./scripts/alignment-filter ../../data/phylogenetic/gisaid-$DATE_SEQ-unsampled.fa <(grep "$X$" results/gisaid-$DATE_TREE-unsampled-subtree.txt | cut -f 1) > results/gisaid-$DATE_TREE-$X-unseq0.fa
#done

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
./scripts/alignment-filter ../../../data/data/mmsa_20210622_masked.fa <(cat results/mmsa-$DATE_TREE-unsample-ids.txt results/gisaid-$DATE_TREE-all-samples0.txt) > ../../data/phylogenetic/mmsa-$DATE_SEQ-sampled-unsampled0.fa 

for X in $SUB_TREES; do
	#waiting for B.1.1.70 B.1.1.7 B.1.177 
	#done C B.1.1.519 B.1.1.317 B.1.221 B.1.36 B.1.351 B.1.160 A
	#for X in ; do
	echo "Enriching tree $X"
	cp results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus
	./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation false --single-child false --ilabel true --set-tip-location false
	./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE-rich.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation true --single-child false --ilabel true --set-tip-location false

	#xzcat results/fasttree/$DATE_TREE-$X/gisaid-$DATE_TREE-ft.xz > results/beast/unsampled/alignment-$DATE_TREE-$X.al
	#mkdir -p results/beast/unsampled/alignment-$DATE_TREE-$X-epa/
	#epa-ng --split results/beast/unsampled/alignment-$DATE_TREE-$X.aln results/beast/unsampled/gisaid-$DATE_TREE-$X-unseq-seq.aln --out-dir results/beast/unsampled/alignment-$DATE_TREE-$X-epa/ --redo
	##fasta uses JC with CAT/Gamma20 (the option used by tree_ft of sarscov2 pipeline
	#
	#epa-ng --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --ref-msa results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference.fasta --query results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query.fasta --out-dir results/beast/unsampled/alignment-$DATE_TREE-$X-epa/ --model JC+G20 --redo
	#
	#scripts/gappa-examine-graft --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --jplace results/beast/unsampled/alignment-$DATE_TREE-$X-epa/epa_result.jplace --out results/beast/run/tree-rich/$X-DTA-$DATE_TREE.MCC.tree -l \"$STATE\" \"non$STATE\" --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv results/gisaid-$DATE_TREE-metadata-unsampled.tsv
	#
	##result: results/beast/run/tree-rich/$X-DTA-$DATE_TREE.MCC.tree


	./scripts/alignment-filter ../../data/phylogenetic/mmsa-$DATE_SEQ-sampled-unsampled0.fa <( grep "$X$" results/gisaid-$DATE_TREE-unsampled-subtree.txt | cut -f 1 ) > results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa 
	./scripts/alignment-filter ../../data/phylogenetic/mmsa-$DATE_SEQ-sampled-unsampled0.fa results/gisaid-$DATE_TREE-$X-samples0.txt > results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa 
	epa-ng --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --ref-msa results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa --query results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa --out-dir results/beast/unsampled/alignment-$DATE_TREE-$X-epa/ --model JC+G20 --redo --baseball-heur --verbose
	scripts/gappa-examine-graft --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --jplace results/beast/unsampled/alignment-$DATE_TREE-$X-epa/epa_result.jplace --out results/beast/run/tree-rich/$X-DTA-$DATE_TREE.MCC.tree -l \"$STATE\" \"non$STATE\" --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv results/gisaid-$DATE_TREE-metadata-unsampled.tsv

done

# for B.1.1.7 extract lineages ourself and put them on files (should be tested, once written and tested, it is deleted and written again)
X=B.1.1.7
mv results/beast/run/tree-rich/$X-DTA-$DATE_TREE.MCC.tree results/beast/run/tree-rich/$X-DTA-$DATE_TREE.MCC.nexus
scripts/split_tree_general --in results/beast/run/tree-rich/$X-DTA-$DATE_TREE.MCC.nexus --out 'results/beast/run/tree-rich/$X-?-DTA-$DATE_TREE.MCC.tree' -l \"$STATE\" \"non$STATE\" --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv 

cd reports-onetree
R
rmarkdown::render('extractLineagesMCC_Germany.Rmd');
rmarkdown::render('extractLineages_State.Rmd');

####
# END: The onetree model: 
####


# on local computer for copying lin-rich
#hzi-get /net/viral_genomics/covid-lineage/huge-lineage-dynamics/analyses/phylogenetic/results/beast/run/lin-rich results/beast/run/


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
mkdir results/beast/run/lin-ius/ results/beast/run/lin-ius/out/
mkdir results/beast/run/lin-ius/aln/ results/beast/run/lin-ius/alnq/
#for X in B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.351 A C B.1.617.2; do
for X in $SUB_TREES; do
	echo "Tree $X "
	#calculate distances between sequences
	#Rscript scripts/filter-by-country.R "$STATE" results/gisaid-$DATE_TREE-metadata-unsampled.tsv results/beast/unsampled/unsampled-$X-in.txt
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
	#Calculating distances between unsampled Germany sequences
	#mash dist -v 0.05 -p 120 results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X.txt.gz
	#zcat results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X.txt.gz | scripts/convert-filename-to-id | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-query-query-$DATE_TREE-$X-id.txt.gz
	
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
	#mkdir -p results/beast/run/lin-ius-3/out/
	#scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius-3/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius-3/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius-3/out/$X-clusters_DTA_MCC_singles_0.5.tsv --allow-new-lineage false --dist-threshold 0.0
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


##DEBUG
#	scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist debug/dist.txt -l \"$STATE\" \"non$STATE\" --out-folder ./debug/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile ./debug/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters ./debug/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples ./debug/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single ./debug/$X-clusters_DTA_MCC_singles.tsv
for X in B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.351 A C B.1.617.2; do
echo $X
#	scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles.tsv
	scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles_0.5.tsv
done
#END DEBUG

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
#lin_ius_clusters_combine lin-ius-3
#lin_ius_clusters_combine lin-ius-no


# Calculate number of samples after sampling and filtering
#for X in $SUB_TREES; do
#./scripts/phylogeo_sankoff_general --in results/gisaid-$DATE_TREE-$X-2.tree --metadata results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --save-nexus --out x --samples y 2>/dev/null
#R -q -e "library(dplyr); library(stringr); samples <- read.table('y'); metadata <- read.table('results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv', sep='\\\\t', quote='\\\"', head=TRUE); metadata[,'country'] <- sapply(strsplit(metadata[,'Location'], '/'), function(x) str_trim(x[2])); data <- inner_join(samples, metadata, by=c('V1' = 'Accession.ID')); data %>% mutate(G = if_else(country == 'Germany', 'Germany', 'nonGermany')) %>% select(G) %>% group_by(G) %>% summarise(n=n());" 2>/dev/null
#done > x2
#cat x2 | grep -v '^>' | grep ' Germany' | awk 'BEGIN{s=0;}{s += $3; }END{print(s);}'
#cat x2 | grep -v '^>' | grep ' nonGermany' | awk 'BEGIN{s=0;}{s += $3; }END{print(s);}'
#

#END OF LIN_IUS

## TESTIG SANKOFF to check number of singletons
mkdir results/beast/run/tree-rich-sk
for X in $SUB_TREES; do
	cp results/beast/run/tree-rich/$X-DTA-$DATE_TREE.MCC.tree results/beast/run/tree-rich-sk/$X-DTA-$DATE_TREE.MCC.tree.nexus
	#./scripts/phylogeo_sankoff_general --in results/beast/run/tree-rich-sk/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/run/tree-rich-sk/$X-DTA-$DATE_TREE.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --print-allow-annotation true --set-tip-location false --save-nexus 
./scripts/phylogeo_sankoff_general --in results/beast/run/tree-rich-sk/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/run/tree-rich-sk/$X-DTA-$DATE_TREE.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --print-allow-annotation true --remove-internal-node-locations --save-nexus --set-tip-location false
done

### Executing another DTA on rich-trees (did not work)
#
#for X in $SUB_TREES; do 
#	./scripts/fill-template.py --template data/X-DTA2-template.xml --name $X --sample <(grep "$X$" results/gisaid-$DATE_TREE-unsampled-subtree.txt | cut -f 1; cat results/gisaid-$DATE_TREE-$X-samples0.txt) --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --out results/beast/$X-DTA2-$DATE_TREE.xml --loc Germany --date $DATE_TREE 
#done
#
#cd /net/viral_genomics/covid-lineage/huge-lineage-dynamics/analyses/phylogenetic/
#for X in $SUB_TREES; do 
#for i in {1..1}; do 
#if [[ "$X" == "B.1.1.7" ]]; then
#qsub -cwd -N DTA2-$X-$i -l h_vmem=64G,mem_free=20G,s_vmem=20G -pe smp 5 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA2.sh run/$X-$i $DATE_TREE $X
#else
#qsub -cwd -N DTA2-$X-$i -l h_vmem=64G,mem_free=20G,s_vmem=20G -pe smp 2 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA2.sh run/$X-$i $DATE_TREE $X
#fi
#done
#done
#
#

## trying with papara
#xzcat results/fasttree/$DATE_TREE-$X/gisaid-$DATE_TREE-ft.xz > results/beast/unsampled/alignment-reference-$X.aln
#scripts/fasta-to-phylip --input-fasta results/beast/unsampled/alignment-reference-$X.aln --output-phy results/beast/unsampled/alignment-reference-$X.phy
#sed 's/EPI_ISL_/S/g' < results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree | sed 's/)[a-z0-9:]*;/);/' > results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.shortname.tree
#
#
#papara -t results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.shortname.tree -s results/beast/unsampled/alignment-reference-$X.phy -q results/gisaid-$DATE_TREE-$X-unseq0.fa -r -n pashmak -j 50
#
#
## trying with hmmer
#hmmbuild2 -g --fast results/beast/unsampled/$X-DTA-$DATE_TREE.0.hmm results/beast/unsampled/alignment-reference-$X.aln
#hmmalign2 --outformat A2M -o results/beast/unsampled/$X-DTA-$DATE_TREE.hmm.aln results/beast/unsampled/$X-DTA-$DATE_TREE.0.hmm results/gisaid-$DATE_TREE-$X-unseq0.fa
#
#epa-ng --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --ref-msa results/beast/unsampled/alignment-reference-$X.aln --query results/beast/unsampled/$X-DTA-$DATE_TREE.hmm.aln --out-dir results/beast/unsampled/alignment-$DATE_TREE-$X-epa/ --model JC+G20 --redo
#
##hmmalign2 --outformat A2M -o test-$X-DTA-$DATE_TREE.hmm.aln results/beast/unsampled/$X-DTA-$DATE_TREE.0.hmm x
##epa-ng --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --ref-msa results/beast/unsampled/alignment-reference-$X.aln --query x --out-dir test-x/ --model JC+G20 --redo
#
#split -l 1000 -d results/gisaid-$DATE_TREE-$X-unseq0.fa results/beast/unsampled/seq-split/$X/gisaid-$DATE_TREE-$X-unseq0.fa.
#
#all:100-287
#triton:
#seq 100 178 | parallel -j 26 -I% --max-args 1 --verbose hmmalign2 --outformat A2M -o results/beast/unsampled/seq-split/$X-aln/$X-DTA-$DATE_TREE.hmm.%.aln results/beast/unsampled/$X-DTA-$DATE_TREE.0.hmm results/beast/unsampled/seq-split/$X/gisaid-$DATE_TREE-$X-unseq0.fa.%
#deimos:
#seq 280 287 | parallel -j 4 -I% --max-args 1 --verbose hmmalign2 --outformat A2M -o results/beast/unsampled/seq-split/$X-aln/$X-DTA-$DATE_TREE.hmm.%.aln results/beast/unsampled/$X-DTA-$DATE_TREE.0.hmm results/beast/unsampled/seq-split/$X/gisaid-$DATE_TREE-$X-unseq0.fa.%
#
#X=B.1.1.7
#for i in {101..287}; do
#qsub -cwd -N ha-$X-$i -l h_vmem=1G,mem_free=1G,s_vmem=1G -o results/beast/unsampled/log/$X-$i -e results/beast/unsampled/log/$X-$i scripts/run-hmmralign2.sh $DATE_TREE $X $i
#done
#
#
##while read a; do
##echo ALIGN HMM: $a;
##done
#
#
### Try mmsa 
#cut -f3 results/gisaid-$DATE_TREE-metadata-unsampled.tsv | tail -n +2  > results/mmsa-$DATE_TREE-unsample-ids.txt 
#
#X=B.1.1.7
#
#cp results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus
#./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation false --single-child false --ilabel true --set-tip-location false
#./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE-rich.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation true --single-child false --ilabel true --set-tip-location false
##
##cut -f3 results/gisaid-$DATE_TREE-metadata-unsampled.tsv | tail -n +2  > results/mmsa-$DATE_TREE-unsample-ids.txt 
##./scripts/alignment-filter ../../../data/data/mmsa_20210622_masked.fa results/mmsa-$DATE_TREE-unsample-ids.txt > ../../data/phylogenetic/mmsa-$DATE_SEQ-unsampled0.fa 
##grep '>' ../../data/phylogenetic/mmsa-$DATE_SEQ-unsampled0.fa | sort | uniq -c | grep -v '^ *1 ' | sort -nr 
###./scripts/rename-alignment-ids results/gisaid-$DATE_TREE-metadata-unsampled.tsv < ../../data/phylogenetic/gisaid-$DATE_SEQ-unsampled0.fa > ../../data/phylogenetic/gisaid-$DATE_SEQ-unsampled.fa
###for X in $SUB_TREES; do
##results/gisaid-$DATE_TREE-$X-samples0.txt > results/gisaid-$DATE_TREE-$X-seq0.fa
#./scripts/alignment-filter ../../../data/data/mmsa_20210622_masked.fa <(cat results/mmsa-$DATE_TREE-unsample-ids.txt results/gisaid-$DATE_TREE-$X-samples0.txt) > ../../data/phylogenetic/mmsa-$DATE_SEQ-sampled-unsampled0.fa 
#./scripts/alignment-filter ../../data/phylogenetic/mmsa-$DATE_SEQ-sampled-unsampled0.fa <(grep "$X$" results/gisaid-$DATE_TREE-unsampled-subtree.txt | cut -f 1) > results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa 
#./scripts/alignment-filter ../../data/phylogenetic/mmsa-$DATE_SEQ-sampled-unsampled0.fa results/gisaid-$DATE_TREE-$X-samples0.txt > results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa 
##	./scripts/alignment-filter ../../data/phylogenetic/mmsa-$DATE_SEQ-unsampled0.fa <(grep "$X$" results/gisaid-$DATE_TREE-unsampled-subtree.txt | cut -f 1) > results/mmsa-$DATE_TREE-$X-unseq0.fa
###done
##epa-ng --split results/beast/unsampled/alignment-$DATE_TREE-$X.aln results/beast/unsampled/gisaid-$DATE_TREE-$X-unseq-seq.aln --out-dir results/beast/unsampled/alignment-$DATE_TREE-$X-epa/ --redo
###fasta uses JC with CAT/Gamma20 (the option used by tree_ft of sarscov2 pipeline
##epa-ng --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --ref-msa results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa --query results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa --out-dir results/beast/unsampled/alignment-$DATE_TREE-$X-epa/ --model JC+G20 --redo
#epa-ng --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --ref-msa results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa --query results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa --out-dir results/beast/unsampled/alignment-$DATE_TREE-$X-epa/ --model JC+G20 --redo --baseball-heur --verbose
## also one is running on ariel
#
###compress (useless)
##scripts/compress-alignment < results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa > results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa-compressed.fa
#
#
######   ADDIGN UNSAMPLD to mcmc trees
