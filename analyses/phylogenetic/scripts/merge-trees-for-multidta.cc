#include "lineage-importation-inject-unsampled.h"
#include <boost/program_options.hpp>

array<string, StateInOut::size> StateInOut::names = {"UK", "nonUK"};

namespace po = boost::program_options;

struct TreeLeafLabelAssignerMergeNexus {
	map<string, int> labels;

    string operator()(string original_label) {
		//cerr << "op()" << original_label << endl;
		if (labels.find(original_label) == labels.end()) {
			labels[original_label] = labels.size() + 1;
		}
		return to_string(labels[original_label]);
		// labels.push_back(original_label);
        // return std::to_string(labels.size());
    }
};


string save_lineage_trees_binary(vector<string>::const_iterator begin, vector<string>::const_iterator end) {
	// assert(begin < lineage_pointer_vector.end());
	// cerr << "save_lineage_trees_binary" << (begin - lineage_pointer_vector.begin()) << " " << (end - lineage_pointer_vector.begin()) << endl;
	assert(begin < end);
	if (end - begin == 1) {
		return *begin;
	}
	return "(" + save_lineage_trees_binary(begin, begin+(end-begin)/2) + "," + 
		save_lineage_trees_binary(begin+(end-begin)/2, end) + "):0.0";
}

void save_lineage_trees(ofstream& fo, TreeLeafLabelAssignerMergeNexus& tree_leaf_label_assigner, const vector<string>& lineage_str) {
	string tree = save_lineage_trees_binary(lineage_str.begin(), lineage_str.end());

	fo << "#NEXUS\n\nBegin taxa;\n\tDimensions ntax=" << tree_leaf_label_assigner.labels.size() << ";\n\tTaxlabels\n";
	vector<pair<string, int>> labels;
	for (auto &l : tree_leaf_label_assigner.labels) {
		labels.push_back(make_pair(l.first, l.second));
	}
	sort(labels.begin(), labels.end(), [&](const pair<string,int> a, const pair<string,int>& b) { return a.second < b.second; });

	for (auto &l : labels) {
		fo << "\t\t" << l.first << endl;
	}
	fo << "\t\t;\nEnd;\n\nBegin trees;\n\tTranslate\n";
	// for (size_t i = 0; i < tree_leaf_label_assigner.labels.size(); i++) {
	bool first = true;
	for (auto &l : labels) {
		if (!first) {
			fo << ",\n";
		}
		first = false;
		fo << "\t\t" << (l.second) << " " << l.first;
	}
	fo << "\n\t\t;\n";
	fo << "tree TREE1 = [&R] " << tree << ";" << endl << "End;" << endl;
}

struct MoveNodesNotInMetadataOutCountry {
	const FullMetadatas& full_metadatas_in_filter;
	MoveNodesNotInMetadataOutCountry(const FullMetadatas& full_metadatas_in_filter) : full_metadatas_in_filter(full_metadatas_in_filter) {
	}

	void run(INode& n) {
		if (n.isLeaf() && full_metadatas_in_filter.data.find(n.label) == full_metadatas_in_filter.data.end()) {
			if (n.data.get_location_in_prob() >= 0.5) {
				cerr << "Not found " << n.label << " " << endl;
			}
			n.data.annot_vars["location.set"] = "{\"" + StateInOut::names[0] + "\",\"" + StateInOut::names[1] + "\"}";
			n.data.annot_vars["location.set.prob"] = "{0.0,1.0}";
		}
		for (INode& c : n.children) {
			run(c);
		}
	}
};


int main(int argc, char* argv[]) {
	// loading nexus files, keep nodes germany with >0.5 probability, merge all remaining trees, move root back in time. save nexus file.

	// Declare the supported options.
	po::options_description desc("Find lineages and inject unsampled to the lineages");
	desc.add_options()
	    ("help", "Run for multiDTA")
	    ("location_label,l", po::value<vector<string>>()->multitoken(), "location labels, e.g. Germany nonGermany")
	    ("trees", po::value<vector<string>>()->multitoken(), "input files")
	    ("root-date", po::value<double>(), "new date for root as a float, e.g. 1900.50")
	    ("metadata", po::value<string>(), "metadata file")
	    ("metadata-in-filter", po::value<string>(), "metadata file, we only consider these as in-country")
	    ("out", po::value<string>()->required(), "output file file")
	    ("log", po::value<string>(), "log file.")
	    ("samples", po::value<string>(), "file for writing samples to.")

	/*
	    ("dist", po::value<string>(), "file of distances")
	    ("treefile", po::value<string>(), "tree file name (only for output)")
	    ("lin-prefix", po::value<string>(), "prefix for lienage names")
	    ("out-folder", po::value<string>(), "output folder")
	    ("out-clusters", po::value<string>(), "output file for clusters")
	    ("out-clusters-single", po::value<string>(), "output file for clusters with single samples")
	    ("out-cluster-samples", po::value<string>(), "output file for cluster samples")
	    ("filter-out", po::value<string>(), "strings to be removed from cluster creation. only for no-tree")
	    ("out-0.5", po::value<vector<string>>()->multitoken(), "cluster, cluster_samples, cluster_Singles files for cutoff=0.5")
	    ("allow-new-lineage", po::value<bool>()->default_value(false), "If this parameter is true, for each query, the most similar sample from STATE or nonSTATE is assigned. So new lineages might be created")
	    ("do-assign", po::value<bool>()->default_value(true), "Assign samples not in the tree based on distance to the inferred lineages.")
	    ("dist-threshold", po::value<double>()->default_value(1), "Distance threshold. All the similarities bellow this value are ignored.")
	    ("no-tree", po::value<bool>()->default_value(false), "If no tree is true, lineages between samples are calculated based on distance file.")
	*/
	;

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
		.options(desc)
		.run();

	po::variables_map vm;
	po::store(parsed_options, vm);

	StateInOut::names = {vm["location_label"].as<vector<string>>()[0], vm["location_label"].as<vector<string>>()[1]};

	cerr << "StateInOut::names " << StateInOut::names << endl;

	string annotation = vm["metadata"].as<string>();
	map<string, Metadata> metadata = load_map(annotation);
	cerr << "metadata loaded " << metadata.size() << endl;

	FullMetadatas full_metadatas;
	full_metadatas.load(vm["metadata"].as<string>());
	FullMetadatas full_metadatas_in_filter;
	full_metadatas_in_filter.load(vm["metadata-in-filter"].as<string>());
	cerr << "metadata-in-filter loaded " << full_metadatas_in_filter.data.size() << " " << vm["metadata-in-filter"].as<string>() << endl;

	double root_date = vm["root-date"].as<double>();


	// output_folder = vm["out-folder"].as<string>();

	////load dist file
	//SampleDistance sample_distance;
	//sample_distance.load(vm["dist"].as<string>());
	//cerr << "sample distance loaded " << sample_distance.dist_.size() << " elements" << endl;

/*
	Lineages lin(metadata);
	double reference_height = -1;

	FullMetadatas full_metadatas;
	full_metadatas.load(vm["metadata"].as<string>());
*/

	ofstream fo_samples;
	if (vm.count("samples") > 0) {
		fo_samples.open(vm["samples"].as<string>());
	}

	ofstream log_file;
	if (vm.count("log") > 0) {
		log_file.open(vm["log"].as<string>());
	}
	
	vector<string> lineage_str;

	TreeLeafLabelAssignerMergeNexus tree_leaf_label_assigner;
	ofstream fo(vm["out"].as<string>());

	for (auto tree_fn : vm["trees"].as<vector<string>>()) {
		Lineages lin(metadata);
		// double reference_height = -1;

		cerr << "Loading tree " << tree_fn << endl;
		INode phylo = load_tree<INode>(tree_fn);

		LineageExtraction alg(metadata, lin);
		alg.assign_date(phylo);

		MoveNodesNotInMetadataOutCountry moveNodesNotInMetadataOutCountry(full_metadatas_in_filter);
		moveNodesNotInMetadataOutCountry.run(phylo);
		
		cerr << " dates assigned to nodes" << endl;
		alg.find_lineages(phylo, 0, "", lin.lineages_half, 0.5, false, log_file);
		cerr << " Lineage found " << lin.lineages_half.size() << endl;
		for (auto& l : lin.lineages_half) {
			if (l.samples.size() == 0) {
				cerr << "Lineage with size=0" << endl;
				continue;
			}
			lineage_str.push_back(lin.lineage_tree(full_metadatas, l, *l.root, 0.5, tree_leaf_label_assigner).first + ":" + to_string(l.root->data.year-root_date));

			if (fo_samples) {
				for (auto &s: l.samples) {
					fo_samples << s << endl;
				}
			}
		}
	}

	save_lineage_trees(fo, tree_leaf_label_assigner, lineage_str);

	cerr << "saved merged " << lineage_str.size() << " lineages to " << vm["out"].as<string>() << endl;

	// INode phylo;

	// if (!vm["no-tree"].as<bool>()) {
	// 	phylo = load_tree<INode>(vm["tree"].as<string>());
	// 	cerr << "tree loaded " << vm["tree"].as<string>() << endl;

	// 	if (vm.count("log") > 0) {
	// 		fill_tree_info(phylo, full_metadatas);
	// 	}

	// 	bool allow_new_lineages = vm["allow-new-lineage"].as<bool>();
	// 	LineageExtraction alg(metadata, lin);
	// 	alg.assign_date(phylo);
	// 	cerr << "dates assigned to nodes" << endl;
	// 	alg.find_lineages(phylo, 0, vm["lin-prefix"].as<string>(), lin.lineages_NA, -1, allow_new_lineages, log_file);
	// 	cerr << "linages found " << lin.lineages_NA.size() << " lienages" << endl;

	// 	alg.find_lineages(phylo, 0, vm["lin-prefix"].as<string>(), lin.lineages_half, 0.5, allow_new_lineages, log_file);
	// 	cerr << "linages found " << lin.lineages_half.size() << " lienages with cutoff=0.5" << endl;

	// 	reference_height = alg.reference_height;

	// 	QueryAssignment queryAssignment(lin);
	// 	if (vm["do-assign"].as<bool>()) {
	// 		//vector<string> query_samples = sample_distance.get_query_samples();
	// 		//alg.assign_to_lineages(query_samples, sample_distance); //, compressing = true
	// 		queryAssignment.assign_to_lineages(vm["dist"].as<string>(), allow_new_lineages, vm["dist-threshold"].as<double>(), metadata, log_file);
	// 		//cerr << "queries assigned to lineages queries=" << query_samples.size() << endl;

	// 		//for (auto& l : alg.lineages) {
	// 		//	cerr << "L " << l.name << " sample:" << l.samples.size() << " query=" << l.query_samples.size() << " tmrca=" << l.tmrca << endl;
	// 		//}

	// 		queryAssignment.remove_empty_lineages();
	// 	}

	// } else {
	// 	QueryClustering queryClustering(lin, vm["lin-prefix"].as<string>());

	// 	if (vm.count("filter-out")) {
	// 		queryClustering.load_filter_out(vm["filter-out"].as<string>());
	// 	}

	// 	queryClustering.assign_to_lineages(vm["dist"].as<string>(), vm["dist-threshold"].as<double>(), metadata);
	// }

	// //save ... 
	// lin.save_clusters(vm["out-clusters"].as<string>(), reference_height, vm["treefile"].as<string>(), vm["out-clusters-single"].as<string>(), lin.lineages_NA);
	// lin.save_cluster_samples(vm["out-cluster-samples"].as<string>(), full_metadatas, lin.lineages_NA);

	// if (vm.count("out-0.5")) {
	// 	vector<string> out_args = vm["out-0.5"].as<vector<string>>();
	// 	string f_clusters = out_args[0], f_cluster_samples = out_args[1], f_clusters_single = out_args[2];
	// 	lin.save_clusters(f_clusters, reference_height, vm["treefile"].as<string>(), f_clusters_single, lin.lineages_half, "0.5");
	// 	lin.save_cluster_samples(f_cluster_samples, full_metadatas, lin.lineages_half);
	// 	cerr << "cutoff lienages saved " << f_clusters << " singles=" << f_clusters_single << " samples=" << f_cluster_samples << endl;
	// 	if (log_file) {
	// 		lin.save_lineage_trees(log_file, full_metadatas, lin.lineages_half, 0.5);

	// 		log_file << "Trees without unsampled" << endl;
	// 		deque<map<string, vector<string>>> lineages_query_of_samples;
	// 		for (auto& l : lin.lineages_half) {
	// 			lineages_query_of_samples.push_back(l.query_of_samples);
	// 			l.query_of_samples = map<string, vector<string>>();
	// 		}
	// 		lin.save_lineage_trees(log_file, full_metadatas, lin.lineages_half, 0.5);
	// 		for (auto& l : lin.lineages_half) {
	// 			l.query_of_samples = lineages_query_of_samples.front();
	// 			lineages_query_of_samples.pop_front();
	// 		}
	// 	}
	// }

	return 0;
}
