# Importation lineage inference 
The computational workflow for inference of importation lineages of SARS-Cov-2 into Germany. The workflow was adapted from [https://github.com/COG-UK/UK-lineage-dynamics-analysis], for processing a 20-fold larger dataset (original workflow for analysis of first UK wave processed 50K samples and this analysis was handling 1.8M sequence samples). 

Goliaei S, Foroughmand-Araabi MH, Roddy A, McHardy A C, Weber A, Oeversti S, Kuehnert D

Note that because of the GISAID terms of use genomic sequences cannot be shared in this repository. Instead, in the metadata files, we kept only ID of sequences from GISAID data.


# Workflow overview:
File `analyses/phylogenetic/run-2.sh` [run-2.sh](analyses/phylogenetic/run-2.sh) is the script containing all the steps. 
The workflow includes some optional steps, some steps for testing the results, and some codes to be executed on HZI servers. Some steps are very time-consuming, so one may want to reproduce other steps. 

## Setting up the environment:
### Required datasets.
We assume that following data are located at:
* Raw sequences at `../../../data/data/gisaid-$DATE_SEQ-raw.fa`.
* Metadata file at `../../../data/data/gisaid-$DATE_METADATA-metadata.tsv`.
* Alignment of sequences at `../../../data/data/mmsa_20210622_masked.fa`

Note that the paths are relative to the working directory [`analyses/phylogenetic/`](analyses/phylogenetic/).

### Required applications
Install the following applications:
* Thorney BEAST
* treeannotator from the BEAST package.
* mash application to be installed in the PATH.
* `boost_program_options` and `boost_iostreams` located on the `$(CONDA_PREFIX)/lib/` folder for building executables.

Note: Environment variable `PATH` should contain above mentioned installed applications.
We highly recommend installation of conda and creation of an environment for the current repository containing the installed applications.

### Building executables from the current repo.
Applications from the current repository should be compiled and executable files produces. Code of some applications are not located in the current repository. So we assume that a clone of repository `https://github.com/hzi-bifo/phylogeo-tools` is located at `../../../phylogeo-tools` (relative to the working directory). Then you can check if the symlinked files are pointing to the current source codes (e.g. `scripts/state.h` file relative to the working directory).

To create executables, on the folder [`scripts/`](analyses/phylogenetic/scripts/) it is enough to run following command:
```
make
```
Note that the [makefile](analyses/phylogenetic/scripts/makefile) uses environment variable `CONDA_PREFIX` based on which include directory `$(CONDA_PREFIX)/include` and lib directory `$(CONDA_PREFIX)/lib` are addressed. 



    
## Setting environment variables:
* Environment variable `LD_LIBRARY_PATH` should point to a place with libboost_programoptions.
* All the executions could be done from working directory [`analyses/phylogenetic/`](analyses/phylogenetic/) relative to the home of the github repository.
* Following variables are reporesenting current data timestamps. So, please keep it as so.
```
DATE_TREE=20210602
DATE_METADATA=20210602
DATE_SEQ=20210602
```
* State should be "Germany". For other countries, some modifications should be made.
```
STATE=Germany
```
* `TWD` points to a place with temporary but large free space.
```
TWD=/net/sgi/viral_genomics/hadi/tmp/
```
* `SUB_TREES` are the important subtrees found after subsampling.
```
SUB_TREES="A B.1.1.7 B.1.1.519 B.1.1.70 B.1.1.317 B.1.177 B.1.160 B.1.221 B.1.36 B.1.258 B.1.351 C"
```

## Cleaning:
* Input (not in repository): `../../../data/data/gisaid-$DATE_SEQ-raw.fa`
* Output (not in repository): `../../../data/data/gisaid-$DATE_SEQ-Wu1.fa` `$TWD/gisaid-$DATE_SEQ/gisaid-$DATE_SEQ-raw.fa` `../../../data/data/gisaid-$DATE_SEQ.fa`
Note that the file `../../../data/data/gisaid-$DATE_SEQ-raw.fa` should be edited to contain only sequence `hCoV-19/Wuhan/WH04/2020` after execution of the first step.
```
#cleaning
grep -A1000 -F 'hCoV-19/Wuhan/WH04/2020' ../../../data/data/gisaid-$DATE_SEQ-raw.fa > ../../../data/data/gisaid-$DATE_SEQ-Wu1.fa
#edit then ..., label: hCoV-19/Wuhan/WH04/2020|EPI_ISL_406801|2020-01-05
DIR=`pwd`
mkdir $TWD/gisaid-$DATE_SEQ/
cd $TWD/gisaid-$DATE_SEQ/
ln -s $DIR/../../../data/data/gisaid-$DATE_SEQ-raw.fa .
bash $DIR/sarscov2phylo/clean_gisaid.sh -i gisaid-$DATE_SEQ-raw.fa -o gisaid-$DATE_SEQ.fa -t 30
mv gisaid-$DATE_SEQ.fa $DIR/../../../data/data/

cd $DIR
```

## Subsampling
* Input: `../../../data/data/gisaid-$DATE_METADATA-metadata.tsv` `../../../data/data/gisaid-$DATE_SEQ.fa`
* Output: `results/pangolin-$DATE_TREE-r0.tree` `results/gisaid-$DATE_TREE-metadata-sampled.tsv` `../../data/phylogenetic/gisaid-$DATE_SEQ-sampled0.fa` `../../data/phylogenetic/gisaid-$DATE_SEQ-sampled.fa` `results/gisaid-$DATE_TREE-?-samples0.txt` (results/gisaid-20210602-B.1.1.7-samples0.txt, ...) (for different subtrees with their names substitutes of `$SUB_TREES` character '?') `results/gisaid-$DATE_TREE-?-sampled0.tree` (results/gisaid-$DATE_TREE-B.1.1.7-sampled0.tree, ...) (for different subtrees of `$SUB_TREES` with their names substitutes character '?') 
```
## SUBSAMPLING: 
#V3 (CURRENT): build tree as pangolin tree: 
./scripts/tree-build-as-pangolin --out results/pangolin-$DATE_TREE-r0.tree --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv --seqs ../../../data/data/gisaid-$DATE_SEQ.fa  -l Germany nonGermany
# we use following for filtering non-good samaple (e.g. sample without complete date)
./scripts/phylogeo_sankoff_general --in results/pangolin-$DATE_TREE-r0.tree --out results/gisaid-$DATE_TREE-$STATE.tree --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv --location_label $STATE non$STATE --cond 2 "==" Germany

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

# Extracting data of the SUB_TREES
./scripts/partition-by-name --in results/gisaid-$DATE_TREE-sampled.tree --par $SUB_TREES --samples "results/gisaid-$DATE_TREE-?-samples0.txt" --trees "results/gisaid-$DATE_TREE-?-sampled0.tree" -l $STATE non$STATE --print-annotation false --print-internal-node-label false

```

## Creating phylogeny for each sub-tree of SUB_TREES:
* Input: `../../data/phylogenetic/gisaid-$DATE_SEQ-sampled.fa` `results/gisaid-$DATE_TREE-$X-samples0.txt` `results/gisaid-$DATE_TREE-$X-seq0.fa` `../../../data/data/gisaid-$DATE_SEQ-Wu1.fa`
* Output:  `results/gisaid-$DATE_TREE-$X-1.tree` 

```
# Creating sequences for each sub-tree
for X in $SUB_TREES; do
        ./scripts/alignment-filter ../../data/phylogenetic/gisaid-$DATE_SEQ-sampled.fa results/gisaid-$DATE_TREE-$X-samples0.txt > results/gisaid-$DATE_TREE-$X-seq0.fa
        cat ../../../data/data/gisaid-$DATE_SEQ-Wu1.fa >> results/gisaid-$DATE_TREE-$X-seq0.fa
done
```

Following commands creates the phylogenies via [sarscov2phylo https://github.com/roblanf/sarscov2phylo] pipeline.
Execution of following commands are time consuming.
```
for X in $SUB_TREES; do
        DIR=`pwd` ;
        rm -rf results/fasttree/$DATE_TREE-$X/ ;
        mkdir -p results/fasttree/$DATE_TREE-$X/ ;
        ln -s $DIR/results/gisaid-$DATE_TREE-$X-seq0.fa results/fasttree/$DATE_TREE-$X/gisaid-$DATE_TREE-$X-seq0.fa ;
        cp scripts/qrun-fasttree.sh results/fasttree/$DATE_TREE-$X/qrun.sh ;
        cd results/fasttree/$DATE_TREE-$X/ ;
        # Following command was made for submitting a job. To be checked in this version
        bash qrun.sh $DATE_SEQ $DATE_TREE $X ;
        cd - ;
done
```



## MCMC for phylogeny topology reconstruction:
* Input: `results/fasttree/$DATE_TREE-$X/ft_FBP.tree` `results/gisaid-$DATE_TREE-$X-samples0.txt` `results/gisaid-$DATE_TREE-metadata-sampled.tsv` `data/X.fixedRootPrior.skygrid-template-thorney.xml`
* Output: `results/beast/run/$X-$i/` for `$X` in `$SUB_TREES` and `$i` in 31-35, `results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.log`
* Intermediate files: `results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.xml` 

Declaring initial clock rate for MCMCs.
```
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
```

Creating XMLs for BEAST.
```
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
```

Executing BEAST MCMCs.
```
for X in $SUB_TREES; do
    #current jobs number are 31..35, previously it was 1..5 for all and 6..20 for B.1.1.7
    for i in {31..35}; do
        # For job sbumission: qsub -cwd -N beast-$X-$i -M foroughmand@gmail.com -l h_vmem=10G,mem_free=10G,s_vmem=10G -pe smp 3 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast.sh run/$X-$i $DATE_TREE $X
        # To be checked:
        bash run-beast.sh run/$X-$i $DATE_TREE $X
    done
done
```

Analysing the logs:
```
for X in $SUB_TREES; do
    for i in {31..35}; do
        loganalyser results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.log
    done
done
```

Combining logs:
```
#combining and making a log file for each subtree
for X in $SUB_TREES; do
    logcombiner `
        for i in {31..35}; do
            echo results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-fixed.log
        done` results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.log
done
```

## Location assignment MCMC Step (DTA)
* Input: `results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.trees` `results/gisaid-$DATE_TREE-metadata-sampled.tsv`
* Output: `results/beast/run/$X-$i/` for `$X` in `$SUB_TREES` and `$i` in 31-32. 
* Intermediate files: `results/beast/$X-DTA-$DATE_TREE.xml` `results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.trees`

Resampling from sampled MCMC trees for next MCMC steps:
```
for X in $SUB_TREES; do
    MCMC_FRP_LOG_FILES=""
    RANGE=`seq 31 35`
    for i in $RANGE; do
        python scripts/resample.py --burnin 15000000 --rate 100000 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sampled.trees
        python scripts/resample.py --burnin 15000000 --count 500 --tree results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE.trees --out results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sub500.trees
        MCMC_FRP_LOG_FILES="$MCMC_FRP_LOG_FILES results/beast/run/$X-$i/$X.fixedRootPrior.skygrid-$DATE_TREE-sub500.trees"
    done
    logcombiner -trees $MCMC_FRP_LOG_FILES results/beast/$X.fixedRootPrior.skygrid-$DATE_TREE.trees
done
```

Generating XML files:
```
for X in $SUB_TREES; do
        ./scripts/fill-template.py --template data/X-DTA-template.xml --name $X --sample results/gisaid-$DATE_TREE-$X-samples1.txt --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/beast/$X-DTA-$DATE_TREE.xml --loc Germany --date $DATE_TREE
done
```

Executing BEAST MCMC:
```
for X in $SUB_TREES; do
    for i in {31..32}; do
        # Submission command: qsub -cwd -N DTA-$X-$i -l h_vmem=64G,mem_free=20G,s_vmem=20G -pe smp 5 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-beast-DTA.sh run/$X-$i $DATE_TREE $X
        # To be checked
        bash run/$X-$i $DATE_TREE $X
    done
done
```

## Creating MCC tree:
* Input: `results/beast/run/$X-$i/$X-DTA-$DATE_TREE.sub4500.trees` for `$X` in `$SUB_TREES` and `$i` in 31-32 `results/beast/run/$X-$i/$X-DTA-$DATE_TREE.log`
* Output: `results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree` for `$X` in `$SUB_TREES`

Combining logs of DTA MCMC:
```
for X in $SUB_TREES; do
    MCMC_DTA_LOG_FILES=""
    MCMC_DTA_TREE_FILES=""
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
```

Running tree-annotator from BEAST command-line tools.
```
mkdir results/beast/run/all
for i in {1..1}; do
    for X in $SUB_TREES; do
        # Job submission command: qsub -cwd -N ta-$X-$i -l h_vmem=64G,mem_free=20G,s_vmem=20G -pe smp 2 -o results/beast/run/out/$X-$i -e results/beast/run/out/$X-$i run-treeannotator.sh $DATE_TREE $X $i
        # To be checked
        bash run-treeannotator.sh $DATE_TREE $X $i
    done
done
```

(Optional) Compressing log files and sampled trees.
```
for X in $SUB_TREES; do  
    xz -f -k -v -T0 results/beast/run/all/$X-DTA-$DATE_TREE.combined.trees
    xz -f -v -T0 results/beast/run/all/$X-DTA-$DATE_TREE.combined.log
done

```

## Generating reports for MCC tree for sub-sampled sequences (Optional)
For the sub-sampled sequences not the results are at `../results/beast/run/all/` (relative to the working directory). 
The reports on the [`reports`](analyses/phylogenetic/reports/) folder could be executed and lineages for this tree could be generated and evaluated.

## Extracting importation lineages and adding unsampled sequences from Germany to the importation lienages
* Input: `results/gisaid-$DATE_TREE-samples.txt` `results/gisaid-$DATE_TREE-unsampled.txt` `results/gisaid-$DATE_TREE-unsampled-subtree.txt` `../../../data/data/gisaid-$DATE_METADATA-metadata.tsv` `../../../data/data/mmsa_20210622_masked.fa`
* Output: `results/beast/run/lin-ius/clusters_DTA_MCC_NA.tsv` `results/beast/run/lin-ius/clusterSamples_DTA_MCC_NA.tsv` `results/beast/run/lin-ius/clusters_DTA_MCC_0.5.tsv` `results/beast/run/lin-ius/clusterSamples_DTA_MCC_0.5.tsv`
* Intermediate files: `results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_NA.tsv` `results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_NA.tsv` `results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles.tsv` `results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_0.5.tsv` `results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_0.5.tsv` `results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles_0.5.tsv`

```
scripts/extract-unsampled --location_label $STATE non$STATE --in results/gisaid-$DATE_TREE-$STATE.tree --samples results/gisaid-$DATE_TREE-samples.txt --unsampled results/gisaid-$DATE_TREE-unsampled.txt $SUB_TREES --cond 2 "==" Germany --metadata ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv
cat results/gisaid-$DATE_TREE-unsampled.txt | grep -v "NA$" > results/gisaid-$DATE_TREE-unsampled-subtree.txt
./scripts/extract-metadata.R results/gisaid-$DATE_TREE-unsampled-subtree.txt ../../../data/data/gisaid-$DATE_METADATA-metadata.tsv results/gisaid-$DATE_TREE-metadata-unsampled.tsv
cat results/gisaid-$DATE_TREE-metadata-sampled.tsv results/gisaid-$DATE_TREE-metadata-unsampled.tsv  > results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv
```

Creating alignment file for unsampled sequences from Germany.
```
cut -f3 results/gisaid-$DATE_TREE-metadata-unsampled.tsv | tail -n +2  > results/mmsa-$DATE_TREE-unsample-ids.txt
for X in $SUB_TREES; do cat results/gisaid-$DATE_TREE-$X-samples0.txt; done > results/gisaid-$DATE_TREE-all-samples0.txt
./scripts/alignment-filter ../../../data/data/mmsa_20210622_masked.fa <(cat results/mmsa-$DATE_TREE-unsample-ids.txt results/gisaid-$DATE_TREE-all-samples0.txt) > ../../data/phylogenetic/mmsa-$DATE_SEQ-sampled-unsampled0.fa

for X in $SUB_TREES; do
        echo "Enriching tree $X"
        cp results/beast/run/all/$X-DTA-$DATE_TREE.MCC.tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus
        ./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation false --single-child false --ilabel true --set-tip-location false
        ./scripts/phylogeo_sankoff_general --in results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out results/beast/unsampled/$X-DTA-$DATE_TREE-rich.MCC.tree --metadata results/gisaid-$DATE_TREE-metadata-unsampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany --nosankoff --print-allow-annotation true --single-child false --ilabel true --set-tip-location false

        ./scripts/alignment-filter ../../data/phylogenetic/mmsa-$DATE_SEQ-sampled-unsampled0.fa <( grep "$X$" results/gisaid-$DATE_TREE-unsampled-subtree.txt | cut -f 1 ) > results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa
        ./scripts/alignment-filter ../../data/phylogenetic/mmsa-$DATE_SEQ-sampled-unsampled0.fa results/gisaid-$DATE_TREE-$X-samples0.txt > results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa
done

```


Some initializations:
```
for X in $SUB_TREES; do
        mv results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fasta
        mv results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fasta
        mv results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa
        mv results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa
done
mkdir results/beast/run/lin-ius/ results/beast/run/lin-ius/out/
mkdir results/beast/run/lin-ius/aln/ results/beast/run/lin-ius/alnq/
```

```
for X in $SUB_TREES; do
        echo "Tree $X "
        #calculate distances between sequences
        csplit -q --prefix=results/beast/run/lin-ius/aln/alignment-$DATE_TREE-$X-epa-reference- --suffix-format=%06d.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/reference-mmsa.fa '/^>/' '{*}'
        rm results/beast/run/lin-ius/aln/alignment-$DATE_TREE-$X-epa-reference-000000.fasta
        csplit -q --prefix=results/beast/run/lin-ius/alnq/alignment-$DATE_TREE-$X-epa-query- --suffix-format=%06d.fasta results/beast/unsampled/alignment-$DATE_TREE-$X-epa/query-mmsa.fa '/^>/' '{*}'
        rm results/beast/run/lin-ius/alnq/alignment-$DATE_TREE-$X-epa-query-000000.fasta
        ls -S results/beast/run/lin-ius/aln/ | grep "alignment-$DATE_TREE-$X-epa-reference-" | awk '{print "'results/beast/run/lin-ius/aln/'" $0}' | grep -v '\-000000.fasta' > results/beast/run/lin-ius/out/alignment-files-$X-reference.txt
        ls -S results/beast/run/lin-ius/alnq/ | grep "alignment-$DATE_TREE-$X-epa-query-" | awk '{print "'results/beast/run/lin-ius/alnq/'" $0}' | grep -v '\-000000.fasta' > results/beast/run/lin-ius/out/alignment-files-$X-query.txt
        mash sketch -o results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-reference.fasta results/beast/run/lin-ius/aln/alignment-$DATE_TREE-$X-epa-reference-*.fasta  2>>e
        mash sketch -o results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta -l results/beast/run/lin-ius/out/alignment-files-$X-query.txt 2>>e
        time mash dist -v 0.05 -p 120 results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-reference.fasta.msh results/beast/run/lin-ius/alignment-$DATE_TREE-$X-epa-query.fasta.msh | ~/bin/pigz-2.6/pigz -k -3 -p10 > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt.gz
        zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X.txt.gz | scripts/convert-filename-to-id | gzip > results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz

        sed 's/Germany+nonGermany/nonGermany/g' results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus > results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus
        scripts/lineage-importation-inject-unsampled --tree results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --dist <( zcat results/beast/run/lin-ius/alignment-distance-$DATE_TREE-$X-id.txt.gz ) -l \"$STATE\" \"non$STATE\" --out-folder results/beast/run/lin-ius/ --metadata results/gisaid-$DATE_TREE-metadata-sampled-unsampled.tsv --lin-prefix LIN-Germany-$X-"$DATE_TREE"_DTA_MCC_ --treefile results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.nexus --out-clusters results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_NA.tsv --out-cluster-samples results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_NA.tsv --out-clusters-single results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles.tsv --out-0.5 results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_0.5.tsv results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_0.5.tsv results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_singles_0.5.tsv

done
```

Creating output files (to be used by R scripts for generating the reports):
```
(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_NA.tsv
for X in $SUB_TREES; do
        tail -n +2 results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_NA.tsv 
done ) > results/beast/run/lin-ius/clusters_DTA_MCC_NA.tsv

(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_NA.tsv
for X in $SUB_TREES; do
        tail -n +2 results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_NA.tsv
done ) > results/beast/run/lin-ius/clusterSamples_DTA_MCC_NA.tsv

(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_0.5.tsv
for X in $SUB_TREES; do
        tail -n +2 results/beast/run/lin-ius/out/$X-clusters_DTA_MCC_0.5.tsv 
done ) > results/beast/run/lin-ius/clusters_DTA_MCC_0.5.tsv

(
X=`echo $SUB_TREES | head -n1 | awk '{print $1;}'`
head -n 1 results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_0.5.tsv
for X in $SUB_TREES; do
        tail -n +2 results/beast/run/lin-ius/out/$X-clusterSamples_DTA_MCC_0.5.tsv
done ) > results/beast/run/lin-ius/clusterSamples_DTA_MCC_0.5.tsv

```

## Internal movement analysis
The internal movement analysis is a Bayesian DTA method that assigns location, state of Germany, to internal nodes of each imported lineage. This is done on the main folders for each subsampling method, e.g. [`analyses/phylogenetic-test-snake-2/`](analyses/phylogenetic-test-snake-2/), [`analyses/phylogenetic-test-subsampling-3/`](analyses/phylogenetic-test-subsampling-3/), and [`analyses/phylogenetic-test-subsampling-5/`](analyses/phylogenetic-test-subsampling-5/).

First, we remove nodes with unknown state.
```
for X in $SUB_TREES; do
  ./scripts/metadata-remove-na.R results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv results/dtamulti/gisaid-$DATE_TREE-$X-metadata-sampled-no-na.tsv
  ./scripts/phylogeo_sankoff_general --in  results/beast/unsampled/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --save-nexus --out results/dtamulti/$X-DTA-$DATE_TREE.MCC.tree.2.nexus --metadata results/gisaid-$DATE_TREE-$X-metadata-sampled.tsv --location_label \"$STATE\" \"non$STATE\" --cond 2 "==" Germany  --print-annotation false --print-internal-node-label false --single-child false --nosankoff --remove-invalid-children true
done

for X in $SUB_TREES; do
  cat results/dtamulti/gisaid-$DATE_TREE-$X-metadata-sampled-no-na.tsv
done > results/dtamulti/gisaid-$DATE_TREE-all-metadata-sampled-no-na.tsv
```

Then, we merge all the imported lineages to one tree with a root with long branches. This helps us in having one set of parameters for all the different lienages.
```
./scripts/merge-trees-for-multidta --trees `for X in $SUB_TREES; do echo results/dtamulti/$X-DTA-$DATE_TREE.MCC.tree.2.nexus; done` --root-date 1900 -l \"Germany\" \"nonGermany\" --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --metadata-in-filter results/dtamulti/gisaid-$DATE_TREE-all-metadata-sampled-no-na.tsv --samples results/dtamulti/sampled.fixedRootPrior.skygrid-$DATE_TREE.samples --out results/dtamulti/sampled.fixedRootPrior.skygrid-$DATE_TREE.trees
```

Then, we fill the template and execute the MCMC via BEAST.
```
./scripts/fill-template.py --template data/X-DTAmulti-template.xml --name sampled --sample results/dtamulti/sampled.fixedRootPrior.skygrid-$DATE_TREE.samples --metadata results/gisaid-$DATE_TREE-metadata-sampled.tsv --out results/dtamulti/sampled-DTAmulti-$DATE_TREE.0.xml --loc Germany --date $DATE_TREE --loc_depth

sed 's/Baden-Württemberg/Baden-Wurttemberg/g' < results/dtamulti/sampled-DTAmulti-$DATE_TREE.0.xml > results/dtamulti/sampled-DTAmulti-$DATE_TREE.xml

for i in {31..32}; do
  bash run-beast-DTAmulti.sh run/sampled-$i $DATE_TREE sampled
done
```

Then, we merge the trees from two independent runs and subsample them for tree annotator.
```
mkdir results/dtamulti/run/all results/dtamulti/run/lin-ius-3/
MCMC_DTA_LOG_FILES=""
MCMC_DTA_TREE_FILES=""
for i in {31..32}; do
srun python scripts/resample.py --burnin 500000 --rate 4500 --tree results/dtamulti/run/sampled-$i/sampled-DTA-$DATE_TREE.trees --out results/dtamulti/run/sampled-$i/sampled-DTA-$DATE_TREE.sub4500.trees
MCMC_DTA_LOG_FILES="$MCMC_DTA_LOG_FILES results/dtamulti/run/sampled-$i/sampled-DTA-$DATE_TREE.log"
MCMC_DTA_TREE_FILES="$MCMC_DTA_TREE_FILES results/dtamulti/run/sampled-$i/sampled-DTA-$DATE_TREE.sub4500.trees"
done
srun logcombiner -trees $MCMC_DTA_TREE_FILES results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.trees
srun logcombiner $MCMC_DTA_LOG_FILES results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.log

loganalyser results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.log > results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.log.analyzed
```

Running tree annotator that finds MCC tree.
```
srun --mem=100G --cpus-per-task=2 --qos verylong --time=24:00:00 /net/viral_genomics/covid-lineage/germany-lineage-dynamics/bin/beast/bin/treeannotator -type mcc results/dtamulti/run/all/sampled-DTA-$DATE_TREE.combined.trees  results/dtamulti/run/all/sampled-DTA-$DATE_TREE.MCC.nexus
```

Finally, the lineages will be extracted for the R report generation scripts
```
srun ./scripts/lineage-importation-extract-multidta --tree results/dtamulti/run/all/sampled-DTA-$DATE_TREE.MCC.nexus --lineage-samples results/beast/run/lin-ius-3/clusterSamples_DTA_MCC_0.5.tsv --out results/dtamulti/run/lin-ius-3/clusterMovement_DTA_MCC_0.5.0.tsv -l \"Baden-Wurttemberg\" \"Bavaria\" \"Berlin\" \"Brandenburg\" \"Bremen\" \"Hamburg\" \"Hesse\" \"Lower Saxony\" \"Mecklenburg-Western Pomerania\" \"North Rhine-Westphalia\" \"Rhineland-Palatinate\" \"Saarland\" \"Saxony\" \"Saxony-Anhalt\" \"Schleswig-Holstein\" \"Thuringia\" --root-date 1900
sed 's/"//g' < results/dtamulti/run/lin-ius-3/clusterMovement_DTA_MCC_0.5.0.tsv > results/dtamulti/run/lin-ius-3/clusterMovement_DTA_MCC_0.5.tsv
```


## Generating the reports
* Input: `results/beast/run/lin-ius/clusters_DTA_MCC_NA.tsv` `results/beast/run/lin-ius/clusterSamples_DTA_MCC_NA.tsv` `results/beast/run/lin-ius/clusters_DTA_MCC_0.5.tsv` `results/beast/run/lin-ius/clusterSamples_DTA_MCC_0.5.tsv`
* Output: `reports-ius/importationSummaryMultiState.pdf` `reports-ius/lineageBreakdownMultiState.pdf` `reports-ius/lineageSummaryMultiState.pdf`

The R reports (Rmd files) on the folder [`reports-ius/`](analyses/phylogenetic/reports-ius/) (relative to working directory) could be executed then. 
The reports use the files generated from the last step (importation lineage extraction and adding unsampled sequences).

## Robustness analysis
We checked robustness of results under different sub-sampling methods.
* Fixed in/out ratio: In this method, a fixed bucket size for each week for Germany is assumed as well as in- and out- ratios. For each week, we performed following steps iteratively. First we picked in-ratio sequences from Germany and out-ratio sequences from non-Germany sequences. If any of these two sets are completely sampled, or number of sampled sequences from Germany reached to the fixed bucket size, the loop is ended. Results of the analysis based on this subsampling are available in [`analyses/phylogenetic-test-subsampling-3/`](analyses/phylogenetic-test-subsampling-3/) for Germany bucket size of 100 and in-ratio and out-ratio both equal to 1. 
  * [`analyses/phylogenetic-test-subsampling-3-nounsampled/`](analyses/phylogenetic-test-subsampling-3-nounsampled/) without identicals
  * [`analyses/phylogenetic-test-subsampling-3-ridentical/`](analyses/phylogenetic-test-subsampling-3-ridentical/) with injected identicals considering non-'n' locus
* Fixed country bucket size (100-25): In this method a fixed size for samples from each country is assumed for each week. If sequences of a country for a week was less than the threshould, all the available sequences for that week is sampled. Results are available in [`analyses/phylogenetic-test-subsampling-5/`](analyses/phylogenetic-test-subsampling-5/) for Germany bucket size value of 100 and all non-Germany countries bucket size of 25. 
  * [`analyses/phylogenetic-test-subsampling-5-nounsampled/`](analyses/phylogenetic-test-subsampling-5-nounsampled/) without identicals
  * [`analyses/phylogenetic-test-subsampling-5-ridentical/`](analyses/phylogenetic-test-subsampling-5-ridentical/) with injected identicals considering non-'n' locus
* Sampling proportional to the number of cases (main method): A fixed size for the number of samples for each week is assumed. Then, this value is divided to the countries on each week proportioanl to their number of cases. The proportion for each country is filled with the sequences, or with all the sequences if number of available sequences is less than this value. Results are available in [`analyses/phylogenetic-test-snake-2/`](analyses/phylogenetic-test-snake-2/) for 10000 as the number of sequences for each week.
  * [`analyses/phylogenetic-test-snake-2-nounsampled/`](analyses/phylogenetic-test-snake-2-nounsampled/) without identicals
  * [`analyses/phylogenetic-test-snake-2-ridentical/`](analyses/phylogenetic-test-snake-2-ridentical/) with injected identicals considering non-'n' locus
