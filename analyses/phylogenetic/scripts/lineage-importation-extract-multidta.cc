#include <boost/program_options.hpp>
#include "tree.h"
//#include "state.h"
#include <tr1/unordered_map>
#include <iomanip>
#include <string>
#include <set>

using namespace std;

//array<string, StateInOut::size> StateInOut::names = {"UK", "nonUK"};

class StateGermanyStates {
public:
/*
	enum Type {
		UNKNOWN=-1, IN=0, OUT=1
	};

	Type type;
	StateGermanyStates() = default;
	constexpr StateGermanyStates(Type _type) : type(_type) { }
	constexpr StateGermanyStates(int _type) : type(static_cast<Type>(_type)) { }
	operator  Type() const { return type; }
	constexpr bool operator == (StateGermanyStates error) const { return type == error.type; }
	constexpr bool operator != (StateGermanyStates error) const { return type != error.type; }    
	constexpr bool operator == (Type errorType) const { return type == errorType; }
	constexpr bool operator != (Type errorType) const { return type != errorType; }
	static const int size = 16;
*/
	int type;

	StateGermanyStates(int type=-1) : type(type) {}

	operator std::string() const {
		if (type >= 0 && (size_t) type < names.size()) {
			return names[type];
		}
		return "";
	}

	static vector<string> names;

	static const int size = 16; //useless here!
	static StateGermanyStates unknown; //useless here!
};
vector<string> StateGermanyStates::names;
StateGermanyStates StateGermanyStates::unknown;

ostream& operator<<(ostream& os, StateGermanyStates s) {
	return os << StateGermanyStates::names[s.type];
}

istream& operator>>(istream& is, StateGermanyStates& st) {
	string s;
	is >> s;
	st.type = -1; // error
	for (size_t i=0; i<StateGermanyStates::names.size(); i++) {
		if (s == StateGermanyStates::names[i]) {
			st.type = i;
		}
	}
	assert(st.type != -1);
	return is;
}


namespace po = boost::program_options;


struct Data;
typedef Node<StateGermanyStates, Data> INode;

struct Data {
	double height, year;
	map<string, string> annot_vars;

	int sample_count_in, sample_count_out;

	vector<string> get_annotation_set(string key) {
		if (annot_vars.find(key) == annot_vars.end()) {
			cerr << "annot_vars not found " << key << endl;
		}
		assert(annot_vars.find(key) != annot_vars.end());
		string val = annot_vars[key];
		if (val.size() > 0 && val[0] == '{') {
			assert(val[0] == '{' && val[val.size() - 1] == '}');
			vector<string> spl = split(val.substr(0, val.size()-1).substr(1), ',');
			// return make_pair(to<T>(spl[0].substr(1)), to<T>(spl[1].substr(0, spl[1].size()-1)));
			return spl;
		} else {
			return vector<string>{val};
		}
	}

	string get_state() {
        // cerr << "annot_vars:" << annot_vars << " " << endl;
		assert(annot_vars.find("location.set") != annot_vars.end() && annot_vars.find("location.set.prob") != annot_vars.end());
		vector<string> lset = get_annotation_set("location.set"),
			lprob = get_annotation_set("location.set.prob");

		double max_p = -1;
		string max_loc = "";
		for (size_t i=0; i<lprob.size(); i++) {
			double p = stod(lprob[i]);
			if (p > max_p) {
				max_p = p;
				max_loc = lset[i];
			}
		}
		//cerr << "get_state " << max_loc << " " << max_p << endl;
		assert(max_loc != "");
		assert(max_p >= 1.0/16);
		return max_loc;
		// string lset = annot_vars.find("location.set")->second,
		//        lprob = annot_vars.find("location.set.prob")->second;
		// if (lset.find(',') != string::npos) {
		// 	pair<double, double> prob_ = get_annotation_pair<double>("location.set.prob");
		// 	// cerr << "prob_ " << prob_.first << " " << prob_.second << endl;
		// 	pair<string, string> lset_ = get_annotation_pair<string>("location.set");
		// 	// cerr << "lset__ " << lset_.first << " " << lset_.second << endl;
		// 	double prob[2] = {prob_.first, prob_.second};
		// 	string lset[2] = {lset_.first, lset_.second};
		// 	for (size_t i=0; i<2; i++) {
        //         // cerr << "     " << lset[i] << " " <<  StateInOut::names[0] << endl;
		// 		if (lset[i] == StateInOut::names[0])
		// 			return prob[i];
		// 	}
		// 	return 0;
		// } else {
		// 	if (lset == "{" + StateInOut::names[0] + "}")
		// 		return 1;
		// 	return 0;
		// }
	}

	void set_annotation_vars(string annotation) {
		//cerr << "TMRCA annot " << annotation << endl;
		map<string, string> r;
		for (size_t i=0; i<annotation.size(); ) {
			string id = "", val="";
			while ( i < annotation.size() && annotation[i] != '=' && annotation[i] != ',') {
				id += annotation[i];
				i++;
			}
			if (i < annotation.size() && annotation[i] == '=') {
				i++;
				bool in_brackets = false;
				while (i < annotation.size() && (annotation[i] != ',' || in_brackets) ) {
					if (annotation[i] == '{') in_brackets = true;
					if (annotation[i] == '}') in_brackets = false;
					val += annotation[i];
					i++;
				}
			}
			if (i < annotation.size())
				i++;
			//cerr << "V: " << id << " " << val << endl;
			r[id] = val;
		}
		//return r;
		annot_vars = r;
	}

	Data() : sample_count_in(0), sample_count_out(0) {
	}

};

struct LineageSample {
	static const string column_name_accession, column_name_cluster;
	vector<string> items;

	LineageSample(const vector<string> items = {}) : items(items) {}
};
const string LineageSample::column_name_accession = "Accession.ID";
const string LineageSample::column_name_cluster = "cluster";


struct LineageSamples {
	map<string, vector<LineageSample>> cluster_samples;
	map<string, string> sample_cluster;

	vector<string> column_names = {"Virus.name", "Type", "Accession.ID", "Collection.date", "Location", "Additional.location.information", "Sequence.length", "Host", "Patient.age", "Gender", "Clade", "Pango.lineage", "Pangolin.version", "Variant", "AA.Substitutions", "Submission.date", "Is.reference.", "Is.complete.", "Is.high.coverage.", "Is.low.coverage.", "N.Content", "GC.Content", "sample_date", "decimal_date", "taxon_label", "is_sampled", "cluster"};
	void load(string fn) {
		cluster_samples.clear();
		ifstream fi(fn);
		int line_cnt = 0, column_index_accession = -1, column_index_cluster = -1;
		for (string line; getline(fi, line); line_cnt++) {
			vector<string> line_split = split(line, '\t');
			if (line_cnt == 0) {
				column_names = line_split;
				column_index_accession = indexOf(line_split, {LineageSample::column_name_accession});
				column_index_cluster = indexOf(line_split, {LineageSample::column_name_cluster});
			} else {
				cluster_samples[line_split[column_index_cluster]].push_back(LineageSample(line_split));
				sample_cluster[line_split[column_index_accession]] = line_split[column_index_cluster];
			}
		}
	}
};

struct Movement {
	string from, to, from_label, to_label;
	//for root -> a child of root height is 0
	double from_year, to_year;
	int height;
	bool is_leaf;


	Movement(string from = "", string to = "", string from_label = "", string to_label = "", double from_year=-1, double to_year = -1, int height = -1, bool is_leaf = false) : from(from), to(to), from_label(from_label), to_label(to_label), from_year(from_year), to_year(to_year), height(height), is_leaf(is_leaf) {
	}

	static string header() {
		return "src\tdst\theight\tis_leaf\tsrc_label\tdst_label\tsrc_date\tdst_date";
	}

	string row(char sep = '\t') const {
		return from + sep + to + sep + to_string(height) + sep + to_string(is_leaf) + sep + from_label + sep + to_label + sep + to_string(from_year) + sep + to_string(to_year);
	}
};

ostream& operator<<(ostream& os, const Movement& m) {
	return os << m.row(' ');
}

struct MovementExtractor {

	vector<pair<vector<string>, vector<Movement>>> cluster_movements;

	int internal_node_names = 0;
	void init_tree(INode& n, double time) {
		if (n.label == "") {
			n.label = "i" + to_string(internal_node_names+1);
			internal_node_names++;
		}
		n.data.set_annotation_vars(n.annotation);
		n.data.year = time;
		for (auto& c: n.children) {
			init_tree(c, time + c.branch_length);
		}
	}

	void dfs(INode& n, bool starting_part, int h, vector<string>& current_samples, vector<Movement>& current_movements) {
		cerr << "dfs " << n.label << " " << starting_part << " " << h << current_samples << " " << endl;
		//if (n.label == "EPI_ISL_1645223") {
		//	cerr << "NODE FOUND! " << n.label << " " << starting_part << " " << h << " " << current_samples << " " << current_movements << endl;
		//}
		if (n.isLeaf()) {
			current_samples.push_back(n.label);
		}
		 
		for (auto& c: n.children) {
			if (starting_part) {
				if(c.branch_length > 0) {
					vector<string> samples;
					vector<Movement> movements;
					dfs(c, false, 0, samples, movements);
					//entring to the lienage
					movements.push_back(Movement("", c.data.get_state(), "", c.label, n.data.year, c.data.year, -1, c.isLeaf()));
					cluster_movements.push_back(make_pair(samples, movements));
					//if (c.label == "EPI_ISL_1645223") {
					//	cerr << "NODE PARENT " << n.label << " " << starting_part << " " << h << " " << samples << " " << movements << endl;
					//}
				} else {
					dfs(c, starting_part, h+1, current_samples, current_movements);
				}
			} else {
				//cerr << "ME " << n.data.get_state() << "->" << c.data.get_state() << " " << h << " " << c.isLeaf() << endl;
				current_movements.push_back(Movement(n.data.get_state(), c.data.get_state(), n.label, c.label, n.data.year, c.data.year, h, c.isLeaf()));
				dfs(c, starting_part, h+1, current_samples, current_movements);
			}
		}
	}

	void run(INode& n) {
		// root should always be outside of the lineages.
		cluster_movements.clear();
		vector<string> samples;
		vector<Movement> movements;
		dfs(n, true, -1e5, samples, movements);
	}
};


int main(int argc, char* argv[]) {
	po::options_description desc("After multidta, find movement of lineages inside states.");
	desc.add_options()
	    ("help", "help")
	    ("tree", po::value<string>()->required(), "input files")
	    ("lineage-samples", po::value<string>()->required(), "File of lineage samples.")
	    ("out", po::value<string>()->required(), "output file file")
	    ("log", po::value<string>(), "log file.")
	    ("location_label,l", po::value<vector<string>>()->multitoken(), "location labels, e.g. Germany nonGermany")
	    ("root-date", po::value<double>(), "time of the root as a float variable, like 1902.03")
	/*
	    ("samples", po::value<string>(), "file for writing samples to.")
	    ("metadata", po::value<string>(), "metadata file")
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

	StateGermanyStates::names = vm["location_label"].as<vector<string>>();

	// string annotation = vm["metadata"].as<string>();
	// map<string, Metadata> metadata = load_map(annotation);
	// cerr << "metadata loaded " << metadata.size() << endl;

	// FullMetadatas full_metadatas;
	// full_metadatas.load(vm["metadata"].as<string>());

	// double root_date = vm["root-date"].as<double>();


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

	double root_date = vm["root-date"].as<double>();

	INode phylo = load_tree<INode>(vm["tree"].as<string>());
	LineageSamples lineage_samples;
	lineage_samples.load(vm["lineage-samples"].as<string>());
	MovementExtractor movement_extractor;
	movement_extractor.init_tree(phylo, root_date);
	movement_extractor.run(phylo);

	int all_movements_count = 0;
	ofstream fo(vm["out"].as<string>());
	fo << "cluster" << "\t" << Movement::header() << endl;
	for (auto& sm: movement_extractor.cluster_movements) {
		vector<string>& samples = sm.first;
		vector<Movement>& movements = sm.second;

		string cl = "";
		for (auto& s : samples) {
			if (lineage_samples.sample_cluster.find(s) == lineage_samples.sample_cluster.end()) {
				cerr << "Bad lineage_samples not found " << s << endl;
			}
			assert(lineage_samples.sample_cluster.find(s) != lineage_samples.sample_cluster.end());
			if (cl == "")
				cl = lineage_samples.sample_cluster[s];
			if (cl != lineage_samples.sample_cluster[s]) {
				cerr << "Bad lineage_samples " << cl << " " << lineage_samples.sample_cluster[s] << " " << s << " " << samples << endl;
			}
			assert(cl == lineage_samples.sample_cluster[s]);
			assert(cl != "");
		}
		
		for (auto& m : movements) {
			fo << cl << "\t" << m.row() << endl;
		}

		all_movements_count += movements.size();

		// cerr << "cm " << samples << " " << cl << endl;
	}

	cerr << "Movements for " << movement_extractor.cluster_movements.size() << " clusters and " << all_movements_count << " movements saved to " << vm["out"].as<string>() << endl;

	return 0;
}
