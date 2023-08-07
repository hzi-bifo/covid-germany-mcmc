#include "tree.h"
//#include "state.h"
#include "json.hpp"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using nlohmann::json;


class StateInOut {
public:
	enum Type {
		UNKNOWN=-1, IN=0, OUT=1
	};

	Type type;
	StateInOut() = default;
	constexpr StateInOut(Type _type) : type(_type) { }
	constexpr StateInOut(int _type) : type(static_cast<Type>(_type)) { }
	operator  Type() const { return type; }
	constexpr bool operator == (StateInOut error) const { return type == error.type; }
	constexpr bool operator != (StateInOut error) const { return type != error.type; }    
	constexpr bool operator == (Type errorType) const { return type == errorType; }
	constexpr bool operator != (Type errorType) const { return type != errorType; }
	static const StateInOut unknown;
	static const StateInOut def;
	static const int size = 2;

	//IN, OUT
	static array<string, 2> names;
};

const StateInOut StateInOut::unknown = StateInOut::Type::UNKNOWN;
const StateInOut StateInOut::def = StateInOut::Type::OUT;
ostream& operator<<(ostream& os, StateInOut s) {
	return os << ((s.type == StateInOut::Type::IN) ? StateInOut::names[0] : (s.type == StateInOut::Type::UNKNOWN ? "U" : StateInOut::names[1])) ;
}

istream& operator>>(istream& is, StateInOut& st) {
	string s;
	is >> s;
	if (s == StateInOut::names[0])
		st.type = StateInOut::Type::IN;
	else if (s == StateInOut::names[1])
		st.type = StateInOut::Type::OUT;
	else {
		st.type = StateInOut::Type::OUT;
		cerr << "W Invalid state " << s << " " << StateInOut::names[0] << " " << StateInOut::names[1] << endl;
		//if (s != "unknown")
		//	cerr << "Invalid state " << s << " " << StateInOut::names[0] << " " << StateInOut::names[1] << endl;
		//assert(s == "unknown");
		//st.type = StateInOut::Type::UNKNOWN;
	}
	return is;
}


double years(const string& date) {
	tm date_tm = {0};
	if (!strptime(date.c_str(), "%Y-%m-%d", &date_tm)) {
		cerr << "Invalid format " << date << endl;
		throw runtime_error("Invalid format " + date);
	}
	mktime(&date_tm);  // t is now your desired time_t

	return (double) date_tm.tm_yday / 366.0 + date_tm.tm_year + 1900;
}       

array<string, StateInOut::size> StateInOut::names = {"UK", "nonUK"};

struct Data {
	bool sampled;
	double height;

	Data() : sampled(false) {}
};
typedef Node<StateInOut, Data> INode;

string replace_all(string s, string from, string to) {
	string r = s;
	for (size_t pos = string::npos; (pos = r.find(from)) != string::npos; ) {
		r.replace(pos, from.size(), to);
	}
	return r;
}

struct TreeDualDFS {

	map<int, INode*> edge_label_parent;
	const map<string, Metadata>& metadata_sampled; 
	const map<string, Metadata>& metadata_unsampled;

	double reference_height;

	TreeDualDFS(const map<string, Metadata>& metadata_sampled, const map<string, Metadata>& metadata_unsampled) :
		metadata_sampled(metadata_sampled), metadata_unsampled(metadata_unsampled) {
	}

	void dfs(INode& n1, INode& n2, INode& parent, double h = 0) {
		if (n1.label == "EPI_ISL_753872")
			cerr << "found " << "EPI_ISL_753872" << endl;
		n1.data.height = h;
		if (n1.isLeaf()) {
			assert(metadata_sampled.find(n1.label) != metadata_sampled.end());
			reference_height = years(metadata_sampled.find(n1.label)->second.date) - h;
		}
		//cerr << "D " << n1.label << " " << n2.label << endl;
		edge_label_parent[n2.edge_label] = &parent;
		if (n1.isLeaf() && n1.label != n2.label) {
			cerr << "E L " << n1.label << " " << n2.label << endl;
		}
		assert(!n1.isLeaf() || n1.label == n2.label);
		if (abs(n1.branch_length - n2.branch_length) >= 1e-5 * n1.branch_length && abs(n1.branch_length - n2.branch_length) >= 1e-6) {
			cerr << "E BL " << n1.branch_length << " " << n2.branch_length << endl;
			assert(abs(n1.branch_length - n2.branch_length) < 1e-5 * n1.branch_length || abs(n1.branch_length - n2.branch_length) < 1e-6);
		}
		assert(n1.children.size() == n2.children.size());
		for (auto i1 = n1.children.begin(), i2 = n2.children.begin(); i1 != n1.children.end(); i1++, i2++) {
			dfs(*i1, *i2, n1, h + i1->branch_length);
		}
	}

	string default_annotation = "length_range={BRANCHLENGTH,BRANCHLENGTH},length_95%_HPD={BRANCHLENGTH,BRANCHLENGTH},length=BRANCHLENGTH,location.rate_95%_HPD={BRANCHLENGTH,BRANCHLENGTH},location.rate_median=BRANCHLENGTH,height_median=BRANCHLENGTH,location.set.prob={1.0},height_range={BRANCHLENGTH,BRANCHLENGTH},location.rate=0.0498375841844921,height_95%_HPD={BRANCHLENGTH,BRANCHLENGTH},length_median=BRANCHLENGTH,BRANCHLENGTH0.03709589698173,location=\"Germany\",location.prob=1.0,location.set={\"Germany\"},location.rate_range={0.024498560208070584,0.09036563764378167},height=BRANCHLENGTH";

	void addJPlaces(json& places) {
		for (auto & e : places) {
			//cerr << "placing " << e << endl;
			json connection = e["p"][0];
			//cerr << "  " << connection << endl;
			//cerr << "  " << connection[0].get<int>() << endl;
			if(edge_label_parent.find(connection[0].get<int>()) == edge_label_parent.end()) {
				cerr << "Parent not exists ... " << connection << " " << connection[0].get<int>() << endl;
			}
			assert(edge_label_parent.find(connection[0].get<int>()) != edge_label_parent.end());
			INode* parent = edge_label_parent[connection[0].get<int>()];
			//if (parent->label.size() > 15)
			//	cerr << " PLL" << parent->label << endl;
			//if (e["n"][0] == "EPI_ISL_930860") {
			//	cerr << "  " << e["n"][0].get<string>() << " " << connection[3].get<double>() << " " << parent << endl;
			//	cerr << "  " << parent->label << endl;
			//	cerr << "  " << parent->branch_length << endl;
			//	cerr << "  " << parent->annotation << endl;
			//}
			double branch_length = connection[3].get<double>() + connection[4].get<double>();

			if(metadata_unsampled.find(e["n"][0].get<string>()) == metadata_unsampled.end()) {
				cerr << "sample not found: " << e["n"][0].get<string>() << " " << metadata_unsampled.size() << endl;
			} else {
				//cerr << "sample found: " << e["n"][0].get<string>() << endl;
			}
			assert(metadata_unsampled.find(e["n"][0].get<string>()) != metadata_unsampled.end());
			branch_length = years(metadata_unsampled.find(e["n"][0].get<string>())->second.date) - reference_height - parent->data.height;

			string a = replace_all(default_annotation, "BRANCHLENGTH", to_string(branch_length));
			//TODO: location should be calculated correctly
			INode n(e["n"][0].get<string>(), branch_length, a, list<INode>(), StateInOut::IN, 1, 1, 1);
			parent->children.push_back(n);
			//if (e["n"][0] == "EPI_ISL_930860")
			//	cerr << "  " << parent->children.size() << endl;
			//cerr << "new node " << e["n"][0] << " " << a << endl;
		}
		cerr << "places added " << places.size() << endl;
	}

	void update_stat(INode& n) {
		//TODO:DEBUG
		for (auto & c : n.children)
			update_stat(c);
		n.update_stat();
	}
};

int main(int argc, char* argv[]) {

	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "Run sankoff")
	    ("tree", po::value<string>(), "input tree")
	    ("out", po::value<string>(), "output tree")
	    ("jplace", po::value<string>(), "output tree")
	    ("metadata", po::value<vector<string>>()->multitoken(), "metadata of sampled, metadata of unsampled cases")
	    ("location_label,l", po::value<vector<string>>()->multitoken(), "location labels, e.g. Germany nonGermany")
	;

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
		.options(desc)
		.run();

	po::variables_map vm;
	po::store(parsed_options, vm);

	StateInOut::names = {vm["location_label"].as<vector<string>>()[0], vm["location_label"].as<vector<string>>()[1]};

	INode phylo = load_tree<INode>(vm["tree"].as<string>()) ;


	//load jplace
	ifstream jplace_f(vm["jplace"].as<string>());
	json jplace;
	jplace_f >> jplace;

	string tree_string = jplace["tree"];
	//cerr << "tree_string: " << tree_string.substr(0, 50) << endl;
	istringstream iss(tree_string);
	//cerr << jplace["tree"] << endl;
	INode jplaceTree = load_tree<INode>(iss);

	//cerr << "metadata files " << vm["metadata"].as<vector<string>>()[0] << " " << vm["metadata"].as<vector<string>>()[1] << endl;
	map<string, Metadata> metadata_sampled = load_map(vm["metadata"].as<vector<string>>()[0]);
	map<string, Metadata> metadata_unampled = load_map(vm["metadata"].as<vector<string>>()[1]);

	TreeDualDFS treeDualDfs(metadata_sampled, metadata_unampled);
	treeDualDfs.dfs(phylo, jplaceTree, phylo);

	treeDualDfs.addJPlaces(jplace["placements"]);
	treeDualDfs.update_stat(phylo);



	ofstream fo(vm["out"].as<string>());
	//NodePrinterGeneral<INode> np;
	//np.print(fo, phylo) << ";" << endl;
	cerr << "tree with " << phylo.size << " nodes " << phylo.sample_size << " samples " << phylo.height << " height " << " saved to " << vm["out"].as<string>() << endl;

	save_nexus_tree(fo, phylo);

	return 0;
}
