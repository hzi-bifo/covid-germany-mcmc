#include "tree.h"
#include "state.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

array<string, StateInOut::size> StateInOut::names = {"UK", "nonUK"};


struct Data {
	bool sampled;

	Data() : sampled(false) {}
};
typedef Node<StateInOut, Data> INode;

//it is also in sample-evenly. remove duplicate
struct TreeUnsampledPrinter {
	ofstream& fo;
	const set<string>& special_parents;
	string active_parent;

	bool merge_additional_location;
	std::vector<std::vector<std::string>> conditions;

	const map<string, Metadata>& metadata;

	TreeUnsampledPrinter(ofstream& fo, const set<string>& special_parents, bool merge_additional_location, std::vector<std::vector<std::string>> conditions,
			const map<string, Metadata>& metadata) : 
		fo(fo), special_parents(special_parents), active_parent("NA"),
    		merge_additional_location(merge_additional_location), conditions(conditions),
		metadata(metadata)
	{
	}

	StateInOut::Type check_location_condition(string label, const map<string, Metadata>& metadata) {
		auto i = metadata.find(label);
		vector<string> location = split(i->second.location, '/');
		if (merge_additional_location)
			location.push_back(i->second.location_add);
		//cerr << "checking for location " << location << endl;
		bool in_location = false;
		for (auto & cond : conditions) {
			bool in_loc_cond = false;
			for (size_t i=0; i<cond.size(); i+=3) {
				size_t loc_index = stoi(cond[i])-1;
				if (loc_index >= location.size()) {
					in_loc_cond = false;
					break;
				}
				string l = trim(location[loc_index]);
				if (cond[i+1] == "==" || cond[i+1] == "eq") {
					in_loc_cond = l == cond[i+2];
				} else if (cond[i+1] == ">=" || cond[i+1] == "gt") {
					in_loc_cond = l.find(cond[i+2]) != string::npos;
				}
			}
			//if (trim(location[1]) == "Germany")
			//	cerr << "  cond " << cond << " " << in_loc_cond << endl;
			//if (in_loc_cond) {
			//	cerr << "matched location " << location << " -- " << cond << endl;
			//}
			in_location |= in_loc_cond;
		}
		if (in_location) {
			return StateInOut::Type::IN;
		} else {
			return StateInOut::Type::OUT;
		}
	}

	void visit(const INode& n) {
		if (special_parents.find(n.label) != special_parents.end())
			active_parent = n.label;
		if (n.isLeaf() && !n.data.sampled && metadata.find(n.label) != metadata.end()) {
			if (check_location_condition(n.label, metadata) == StateInOut::IN) 
				fo << n.label << "\t" << active_parent << endl;
		}
	}


	void finish(const INode& n) {
		if (special_parents.find(n.label) != special_parents.end())
			active_parent = "NA";
	}
};

struct TreeSetSampled {
	const set<string>& sampled;

	TreeSetSampled(const set<string>& sampled) : sampled(sampled) {}

	void visit(INode& n) {
		if (sampled.find(n.label) != sampled.end())
			n.data.sampled = true;
	}

	void finish(const INode& n) {
	}
};


int main(int argc, char* argv[]) {

	// Declare the supported options.
	po::options_description desc("Finds leaf nodes of special_subtrees not sampled that satisfy location condition");
	desc.add_options()
	    ("help", "Run sankoff")
	    ("merge", "Merge additional location to the location as its 4rd element")
	    ("metadata", po::value<string>(), "metadata file")
	    ("samples", po::value<string>(), "list of sampled nodes")
	    ("in", po::value<string>(), "input tree")
	    ("out", po::value<string>(), "output tree")
	    ("location_label", po::value<vector<string>>()->multitoken(), "location labels, e.g. Germany nonGermany")
	    ("cost", po::value<vector<float>>()->multitoken()->default_value(vector<float>{0, 1, 1, 0}, "0 1 1 0"), "cost function, in->in, in->out, out->in, out->out")
	    ("cond", po::value<vector<string>>()->multitoken(), "conditions e.g. --cond 2 == Germany 3 >= Dusseldorf --cond 2 == Germany 4 >= Dusseldorf")
	    ("print-annotation", po::value<bool>()->default_value(true), "print annotation")
	    ("print-allow-annotation", po::value<bool>()->default_value(true), "allow printing annotations. If setted to false, annotations are removed.")
	    ("print-internal-node-label", po::value<bool>()->default_value(true), "print internal node labels")
	    ("ilabel", po::value<bool>()->default_value(false), "override internal node labels")
	    ("single-child", po::value<bool>()->default_value(true), "allow single-child internal nodes")
	    //("unsampled", po::value<string>(), "output unsampled file")
	    ("unsampled", po::value<vector<string>>()->multitoken(), "output file containing unsampled samples which are annotated as InState annotated with the first parent from the set. file-name [LIST OF PARENTS TO WHICH SAMPLES ARE ANNOTATED]")
	;

	// Just parse the options without storing them in a map.
	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
		.options(desc)
		.run();
	//cerr << "options run" << endl;

	// Build list of multi-valued option instances. We iterate through
	// each command-line option, whether it is repeated or not. We
	// accumulate the values for our multi-valued option in a
	// container.
	std::vector<std::vector<std::string>> conditions;
	for (const po::option& o : parsed_options.options) {
		if (o.string_key == "cond")
			conditions.push_back(o.value);
	}

	//cerr << "cond created" << endl;

	// If we had other normal options, we would store them in a map
	// here. In this demo program it isn't really necessary because
	// we are only interested in our special multi-valued option.
	po::variables_map vm;
	po::store(parsed_options, vm);

	//cerr << "options stored " << vm["in"].as<string>() << endl;

	StateInOut::names = {vm["location_label"].as<vector<string>>()[0], vm["location_label"].as<vector<string>>()[1]};

	INode phylo = load_tree<INode>(vm["in"].as<string>()) ;

	//cerr << "Location assigned " << in_location_count << " as in " << endl;

	//cerr << "isolate_matrix:" << isolate_matrix << endl;
	//cerr << "cost" << cost << endl;

	//phylo.set_tip_location(isolate_matrix);

	//phylo.annotation= "location=Germany";

	//int removed_count = 0;
	//ofstream fo(splitted_tree_prefix + "1" + ".trees");
	//phylo.remove_invalid_children(id_to_name, removed_count);

	ifstream fi(vm["samples"].as<string>());
	set<string> samples;
	for (string line; fi >> line; ) 
		samples.insert(line);
	cerr << "samples loaded " << samples.size() << " " << *samples.begin() << endl;

	map<string, Metadata> metadata = load_map(vm["metadata"].as<string>());


	TreeSetSampled tss = TreeSetSampled(samples);
	TreeDFSGeneral<INode, TreeSetSampled> setSampled(tss);
	setSampled.dfs(phylo);


	if (vm["ilabel"].as<bool>() == true) {
		InternalNodeLabeler<INode> internalNodeLabeler;
		internalNodeLabeler.run(phylo);
		cerr << "internal nodes relabeled" << endl;
	}
	if (vm["single-child"].as<bool>() == false) {
		SingleChildInternalNodeRemover<INode> singleChildInternalNodeRemover;
		phylo = singleChildInternalNodeRemover.run(phylo);
		cerr << "Single child internal nodes removed " << singleChildInternalNodeRemover.removed_internal_count << endl;
	}


		vector<string> parents = vm["unsampled"].as<vector<string>>();
		parents.erase(parents.begin());
		//sampler.print_unsampled(vm["unsampled"].as<vector<string>>()[0], parents, phylo);
		set<string> special_parents(parents.begin(), parents.end());
		ofstream fo(vm["unsampled"].as<vector<string>>()[0]);
		TreeUnsampledPrinter tup(fo, special_parents, vm.count("merge") > 0, conditions, metadata);
		TreeDFSGeneral<INode, TreeUnsampledPrinter> dfs(tup);

		dfs.dfs(phylo);

	return 0;
}
