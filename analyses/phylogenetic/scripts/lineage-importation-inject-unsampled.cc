#include "lineage-importation-inject-unsampled.h"
#include <boost/program_options.hpp>


namespace po = boost::program_options;

array<string, StateInOut::size> StateInOut::names = {"UK", "nonUK"};

int main(int argc, char* argv[]) {

	// Declare the supported options.
	po::options_description desc("Find lineages and inject unsampled to the lineages");
	desc.add_options()
	    ("help", "Run sankoff")
	    ("tree", po::value<string>(), "input tree")
	    ("dist", po::value<string>(), "file of distances")
	    ("treefile", po::value<string>(), "tree file name (only for output)")
	    ("lin-prefix", po::value<string>(), "prefix for lienage names")
	    ("out-folder", po::value<string>(), "output folder")
	    ("out-clusters", po::value<string>(), "output file for clusters")
	    ("out-clusters-single", po::value<string>(), "output file for clusters with single samples")
	    ("out-cluster-samples", po::value<string>(), "output file for cluster samples")
	    ("metadata", po::value<string>(), "metadata file")
	    ("filter-out", po::value<string>(), "strings to be removed from cluster creation. only for no-tree")
	    ("location_label,l", po::value<vector<string>>()->multitoken(), "location labels, e.g. Germany nonGermany")
	    ("out-0.5", po::value<vector<string>>()->multitoken(), "cluster, cluster_samples, cluster_Singles files for cutoff=0.5")
	    ("allow-new-lineage", po::value<bool>()->default_value(false), "If this parameter is true, for each query, the most similar sample from STATE or nonSTATE is assigned. So new lineages might be created")
	    ("do-assign", po::value<bool>()->default_value(true), "Assign samples not in the tree based on distance to the inferred lineages.")
	    ("dist-threshold", po::value<double>()->default_value(1), "Distance threshold. All the similarities bellow this value are ignored.")
	    ("no-tree", po::value<bool>()->default_value(false), "If no tree is true, lineages between samples are calculated based on distance file.")
	    ("log", po::value<string>(), "log file.")
	;

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
		.options(desc)
		.run();

	po::variables_map vm;
	po::store(parsed_options, vm);

	StateInOut::names = {vm["location_label"].as<vector<string>>()[0], vm["location_label"].as<vector<string>>()[1]};

	string annotation = vm["metadata"].as<string>(),
		output_folder = vm["out-folder"].as<string>();
	map<string, Metadata> metadata = load_map(annotation);
	cerr << "metadata loaded " << metadata.size() << endl;

	////load dist file
	//SampleDistance sample_distance;
	//sample_distance.load(vm["dist"].as<string>());
	//cerr << "sample distance loaded " << sample_distance.dist_.size() << " elements" << endl;

	Lineages lin(metadata);
	double reference_height = -1;

	FullMetadatas full_metadatas;
	full_metadatas.load(vm["metadata"].as<string>());

	ofstream log_file;
	if (vm.count("log") > 0) {
		log_file.open(vm["log"].as<string>());
	}

	INode phylo;

	if (!vm["no-tree"].as<bool>()) {
		phylo = load_tree<INode>(vm["tree"].as<string>());
		cerr << "tree loaded " << vm["tree"].as<string>() << endl;

		if (vm.count("log") > 0) {
			fill_tree_info(phylo, full_metadatas);
		}

		bool allow_new_lineages = vm["allow-new-lineage"].as<bool>();
		LineageExtraction alg(metadata, lin);
		alg.assign_date(phylo);
		cerr << "dates assigned to nodes" << endl;
		alg.find_lineages(phylo, 0, vm["lin-prefix"].as<string>(), lin.lineages_NA, -1, allow_new_lineages, log_file);
		cerr << "linages found " << lin.lineages_NA.size() << " lienages" << endl;

		alg.find_lineages(phylo, 0, vm["lin-prefix"].as<string>(), lin.lineages_half, 0.5, allow_new_lineages, log_file);
		cerr << "linages found " << lin.lineages_half.size() << " lienages with cutoff=0.5" << endl;

		reference_height = alg.reference_height;

		QueryAssignment queryAssignment(lin);
		if (vm["do-assign"].as<bool>()) {
			//vector<string> query_samples = sample_distance.get_query_samples();
			//alg.assign_to_lineages(query_samples, sample_distance); //, compressing = true
			queryAssignment.assign_to_lineages(vm["dist"].as<string>(), allow_new_lineages, vm["dist-threshold"].as<double>(), metadata, log_file);
			//cerr << "queries assigned to lineages queries=" << query_samples.size() << endl;

			//for (auto& l : alg.lineages) {
			//	cerr << "L " << l.name << " sample:" << l.samples.size() << " query=" << l.query_samples.size() << " tmrca=" << l.tmrca << endl;
			//}

			queryAssignment.remove_empty_lineages();
		}

	} else {
		QueryClustering queryClustering(lin, vm["lin-prefix"].as<string>());

		if (vm.count("filter-out")) {
			queryClustering.load_filter_out(vm["filter-out"].as<string>());
		}

		queryClustering.assign_to_lineages(vm["dist"].as<string>(), vm["dist-threshold"].as<double>(), metadata);
	}

	//save ... 
	lin.save_clusters(vm["out-clusters"].as<string>(), reference_height, vm["treefile"].as<string>(), vm["out-clusters-single"].as<string>(), lin.lineages_NA);
	lin.save_cluster_samples(vm["out-cluster-samples"].as<string>(), full_metadatas, lin.lineages_NA);

	if (vm.count("out-0.5")) {
		vector<string> out_args = vm["out-0.5"].as<vector<string>>();
		string f_clusters = out_args[0], f_cluster_samples = out_args[1], f_clusters_single = out_args[2];
		lin.save_clusters(f_clusters, reference_height, vm["treefile"].as<string>(), f_clusters_single, lin.lineages_half, "0.5");
		lin.save_cluster_samples(f_cluster_samples, full_metadatas, lin.lineages_half);
		cerr << "cutoff lienages saved " << f_clusters << " singles=" << f_clusters_single << " samples=" << f_cluster_samples << endl;
		if (log_file) {
			TreeLeafLabelAssignerIdentical treeLeafLabelAssignerIdentical;
			lin.save_lineage_trees(log_file, full_metadatas, lin.lineages_half, 0.5, treeLeafLabelAssignerIdentical);

			log_file << "Trees without unsampled" << endl;
			deque<map<string, vector<string>>> lineages_query_of_samples;
			for (auto& l : lin.lineages_half) {
				lineages_query_of_samples.push_back(l.query_of_samples);
				l.query_of_samples = map<string, vector<string>>();
			}
			lin.save_lineage_trees(log_file, full_metadatas, lin.lineages_half, 0.5, treeLeafLabelAssignerIdentical);
			for (auto& l : lin.lineages_half) {
				l.query_of_samples = lineages_query_of_samples.front();
				lineages_query_of_samples.pop_front();
			}
		}
	}

	return 0;
}
