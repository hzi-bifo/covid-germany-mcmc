ALL: partition sampling test tree-stat alignment-filter partition-by-name phylogeo_sankoff_general sampling-by-name sample-and-partition-by-name sample-evenly tree-build-as-pangolin rename-alignment-ids contract_short_branch extract-unsampled gappa-examine-graft compress-alignment lineage-importation-inject-unsampled convert-filename-to-id filter-alignment-distance tree-stat2 tree-filter-by-metadata partition-by-sample-size merge-trees-for-multidta merge-trees-for-multidta lineage-importation-extract-multidta

partition: partition.cc tree.h
	g++ -o partition partition.cc -O2 -Wall -std=c++11 -lboost_iostreams  -I$(CONDA_PREFIX)/include -L$(CONDA_PREFIX)/lib -Wl,-rpath-link=$(CONDA_PREFIX)/lib -g

sampling: sampling.cc tree.h rangetree.h
	g++ -o sampling sampling.cc -O2 -Wall -std=c++11 -lboost_program_options

test: test.cc rangetree.h
	g++ -o test test.cc -O0 -Wall -std=c++11

tree-stat: tree-stat.cc tree.h
	g++ -o tree-stat tree-stat.cc -O0 -Wall -std=c++11 -lboost_program_options -L$(CONDA_PREFIX)/lib

tree-stat2: tree-stat2.cc tree.h
	g++ -o tree-stat2 tree-stat2.cc -O0 -Wall -std=c++11 -lboost_program_options -L$(CONDA_PREFIX)/lib

partition-by-sample-size: partition-by-sample-size.cc tree.h
	g++ -o partition-by-sample-size partition-by-sample-size.cc -O0 -Wall -std=c++11 -lboost_program_options -L$(CONDA_PREFIX)/lib

tree-filter-by-metadata: tree-filter-by-metadata.cc tree.h state.h
	g++ -o tree-filter-by-metadata tree-filter-by-metadata.cc -O0 -Wall -std=c++11 -lboost_program_options -L$(CONDA_PREFIX)/lib

alignment-filter: alignment-filter.cc 
	g++ -o alignment-filter alignment-filter.cc -O3 -Wall -std=c++11 -lboost_iostreams  -I$(CONDA_PREFIX)/include -L$(CONDA_PREFIX)/lib -Wl,-rpath-link=$(CONDA_PREFIX)/lib

partition-by-name: partition-by-name.cc tree.h
	g++ -o partition-by-name partition-by-name.cc -O0 -Wall -std=c++11 -lboost_iostreams  -I$(CONDA_PREFIX)/include -lboost_program_options -L$(CONDA_PREFIX)/lib

phylogeo_sankoff_general: phylogeo_sankoff_general.cc tree.h
	g++ -o phylogeo_sankoff_general phylogeo_sankoff_general.cc -O2 -Wall -std=c++11 -lboost_iostreams  -I$(CONDA_PREFIX)/include -lboost_program_options

sampling-by-name: sampling-by-name.cc tree.h
	g++ -o sampling-by-name sampling-by-name.cc -O2 -Wall -std=c++11 -lboost_iostreams  -I$(CONDA_PREFIX)/include -lboost_program_options

sample-and-partition-by-name: sample-and-partition-by-name.cc tree.h
	g++ -o sample-and-partition-by-name sample-and-partition-by-name.cc -O2 -Wall -std=c++11 -lboost_iostreams  -I$(CONDA_PREFIX)/include -lboost_program_options

sample-evenly: sample-evenly.cc tree.h rangetree.h
	g++ -o sample-evenly sample-evenly.cc -g -O2 -Wall -std=c++11 -lboost_program_options

tree-build-as-pangolin: tree-build-as-pangolin.cc tree.h 
	g++ -o tree-build-as-pangolin tree-build-as-pangolin.cc -O2 -Wall -std=c++11 -lboost_program_options

rename-alignment-ids: rename-alignment-ids.cc tree.h
	g++ -o rename-alignment-ids rename-alignment-ids.cc -O2 -Wall -std=c++11

contract_short_branch: contract_short_branch.cc tree.h state.h
	g++ contract_short_branch.cc -o $@ -O2 -std=c++11 -Wall -lboost_program_options -L$(CONDA_PREFIX)/lib

extract-unsampled: extract-unsampled.cc tree.h state.h
	g++ extract-unsampled.cc -o $@ -O2 -std=c++11 -Wall -lboost_program_options -g

gappa-examine-graft: gappa-examine-graft.cc tree.h state.h
	g++ gappa-examine-graft.cc -o $@ -O0 -std=c++11 -Wall -lboost_program_options

lineage-importation-inject-unsampled: lineage-importation-inject-unsampled.cc lineage-importation-inject-unsampled.h tree.h state.h
	g++ lineage-importation-inject-unsampled.cc -o $@ -O6 -std=c++11 -Wall -lboost_program_options -I$(CONDA_PREFIX)/include -L$(CONDA_PREFIX)/lib

merge-trees-for-multidta: lineage-importation-inject-unsampled.cc merge-trees-for-multidta.cc lineage-importation-inject-unsampled.h tree.h state.h
	g++ merge-trees-for-multidta.cc -o $@ -O6 -std=c++11 -Wall -lboost_program_options -I$(CONDA_PREFIX)/include -L$(CONDA_PREFIX)/lib

lineage-importation-extract-multidta: lineage-importation-extract-multidta.cc tree.h state.h
	g++ lineage-importation-extract-multidta.cc -o $@ -O6 -g -std=c++11 -Wall -lboost_program_options -I$(CONDA_PREFIX)/include -L$(CONDA_PREFIX)/lib

convert-filename-to-id: convert-filename-to-id.cc 
	g++ convert-filename-to-id.cc -o $@ -O3 -std=c++11 -Wall 

compress-alignment: compress-alignment.cc fasta.h
	g++ compress-alignment.cc -o $@ -O0 -std=c++11 -Wall

filter-alignment-distance: filter-alignment-distance.cc 
	g++ filter-alignment-distance.cc -o $@ -O3 -std=c++11 -Wall 

