#include "tree.h"
#include "state.h"
#include <boost/program_options.hpp>
#include "rangetree.h"
#include <vector>
#include <numeric>
#include <fstream>
#include <set>

namespace po = boost::program_options;

array<string, StateInOut::size> StateInOut::names = {"in", "out"};

struct INodeData {
	int in_count;
	int poc_count, poc_in_count;
	INodeData(int in_count=0, int poc_count=0, int poc_in_count=0): in_count(in_count), poc_count(poc_count), poc_in_count(poc_in_count) {
	}

	INodeData operator+(const INodeData& b) const {
		INodeData r = *this;
		r += b;
		return r;
	}

	INodeData& operator+=(const INodeData& b) {
		in_count += b.in_count;
		poc_count += b.poc_count;
		poc_in_count += b.poc_in_count;
		return *this;
	}
};

typedef Node<StateInOut, INodeData> INode;

vector<string> load_list_file(string fn) {
	vector<string> ret;
	ifstream fi(fn);
	for (string s; std::getline(fi, s); ) {
		ret.push_back(s);
	}
	return ret;
}

struct OutSet {
	vector<string> samples;
	vector<bool> is_sample_poor;

	OutSet(std::initializer_list<string> samples, bool is_poor=false) :samples(samples) {
		assert(samples.size() <= 1);
		for (size_t i=0; i<samples.size(); i++) {
			is_sample_poor.push_back(is_poor);
		}
	}
	OutSet(const std::vector<string>& samples, bool is_poor=false) :samples(samples) {
		for (size_t i=0; i<samples.size(); i++) {
			is_sample_poor.push_back(is_poor);
		}
	}

	//void insert_poors(const OutSet& o) {
	//	samples.insert(samples.end(), o.samples.begin(), o.samples.end());
	//	for (size_t i=0; i<o.samples.size(); i++) {
	//		is_sample_poor.push_back(true);
	//	}
	//}

	void insert_poors(const vector<string>& o) {
		samples.insert(samples.end(), o.begin(), o.end());
		for (size_t i=0; i<o.size(); i++) {
			is_sample_poor.push_back(true);
		}
	}

	operator const vector<string>&() const {
		return samples;
	}

	operator vector<string>&() {
		return samples;
	}

};

ostream& operator<<(ostream& os, const OutSet& o) {
	return os << o.samples;
}

struct Stat {
	int small_size, huge_size, min_in;
	vector<string> selected_nodes;
	set<string> samples_in;

	map<string, StateInOut::Type> sample_location_map;

	Stat(const set<string>& samples_in, int small_size=0, int huge_size = 0, int min_in=0) : small_size(small_size), huge_size(huge_size), min_in(min_in), samples_in(samples_in) {}

	string mul(char c, int m) {
		string r = "";
		for (int i=0; i<m; i++)
			r += c;
		return r;
	}

	// returns outsets
	vector<OutSet> run(INode& n) {
		selected_nodes.clear();
		sample_location_map.clear();
		dfs1(n);
		auto r = dfs(n, 0);
		return get<3>(r);
	}

	//fills node.data
	INodeData dfs1(INode& n) {
		if (n.isLeaf()) {
			sample_location_map[n.label] = n.location;
			return n.data = INodeData(n.location == StateInOut::IN ? 1 : 0, 
				samples_in.find(n.label) != samples_in.end() ? 1 : 0,
				n.location == StateInOut::IN && samples_in.find(n.label) != samples_in.end() ? 1 : 0
				);
		}
		INodeData state_sample_count;
		for (auto & c : n.children)
			state_sample_count += dfs1(c);
		return n.data = state_sample_count;
	}

	vector<string> dfs_all_sampled(const INode& n) {
		vector<string> ret;
		if (n.isLeaf() && n.data.poc_count > 0)
			ret.push_back(n.label);
		for (auto const& c: n.children) {
			vector<string> cr = dfs_all_sampled(c);
			ret.insert(ret.end(), cr.begin(), cr.end());
		}
		return ret;
	}

	int size(const vector<OutSet>& v) const {
		int r = 0;
		for (auto const&vv: v) {
			r += vv.samples.size();
		}
		return r;
	}

	// ret.first = distance to the active set
	// ret.second = active set index
	// ret.third = is the active set set of poors (true) or the closest set (false)
	// ret.fourth = sets
	tuple<int, int, bool, vector<OutSet>> dfs(const INode& n, int d) {
		vector<OutSet> ret;
		if (n.data.poc_count == 0) {
			auto x = dfs_all_sampled(n);
			assert(x.size() == 0);
			ret.push_back(OutSet({}, true));
			return make_tuple(1e5, 0, true, ret);
		}
		if (n.isLeaf()) {
			assert(n.data.poc_count == 1);
			ret.push_back(OutSet({n.label}, min_in > 1 ? true : false));
			return make_tuple(min_in > 1 ? 1e5 : 0, 0, min_in > 1 ? true : false, ret);
		}
		if (n.data.poc_count >= min_in && n.data.poc_count <= huge_size) {
			ret.push_back(OutSet(dfs_all_sampled(n), false));
			return make_tuple(0, 0, false, ret);
		}
		vector<string> poors;
		int min_dist = 1e5, active_index = -1;
		//bool is_active_poor = true;
		for (auto const& c : n.children) {
			auto cr = dfs(c, d+1);
			if (get<2>(cr) == true) {
				assert(get<0>(cr) == 1e5 && get<1>(cr) == 0 && get<3>(cr).size() == 1);
				//ret[active_index].insert(ret[active_index].end(), get<3>(cr)[get<1>(cr)].begin(), get<3>(cr)[get<1>(cr)].end());
				//ret[active_index].insert_poors(get<3>(cr)[get<1>(cr)]);
				// is_active_poor = is_active_poor; if our active is poor, it remains so, otherwise it remains so, too.
				poors.insert(poors.end(), get<3>(cr).back().samples.begin(), get<3>(cr).back().samples.end());
			} else {
				if (min_dist > get<0>(cr)+1) {
					active_index = get<1>(cr) + ret.size();
					min_dist = get<0>(cr) + 1;
				}
				ret.insert(ret.end(), get<3>(cr).begin(), get<3>(cr).end());
				/*
				if (is_active_poor) {
					assert(ret.size() == 1 && active_index == 0 && min_dist == 1e5);
					vector<string> our_poor = ret.back();
					ret = get<3>(cr);
					min_dist = get<0>(cr) + 1;
					active_index = get<1>(cr);
					//ret[active_index].insert(ret[active_index].end(), our_poor.begin(), our_poor.end());
					ret[active_index].insert_poors(our_poor);
					is_active_poor = false;
				} else {
					if (min_dist > get<0>(cr)+1) {
						active_index = get<1>(cr) + ret.size();
						min_dist = get<0>(cr) + 1;
					}
					ret.insert(ret.end(), get<3>(cr).begin(), get<3>(cr).end());
				}*/
			}
		}
		if (ret.size() > 0) {
			//we are not poor
			assert(active_index != -1);
			ret[active_index].insert_poors(poors);
			if (size(ret) != n.data.poc_count) {
				cerr << "E! " << size(ret) << " " << n.data.poc_count << endl;
				for (auto const& c : n.children) {
					auto cr = dfs(c, d+1);
					cerr << " c " << get<3>(cr) << " " << get<0>(cr) << " " << get<1>(cr) << " " << get<2>(cr) << " " << endl;
				}
				cerr << " || poors=" << poors << endl << " || ret=" << ret << endl;
			}
			assert(size(ret) == n.data.poc_count);
			return make_tuple(min_dist, active_index, false, ret);
		} else {
			ret.push_back(OutSet(poors, true));
			assert(size(ret) == n.data.poc_count);
			return make_tuple(min_dist, 0, true, ret);
		}
	}
};

int main(int argc, char* argv[]) {
	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "Run statistics")
	    ("in", po::value<string>(), "input tree")
	    ("in-samples", po::value<string>(), "sampled sequences")
	    ("out-samples", po::value<string>(), "sampled sequences will be placed in files with the given template. Question mark will be replaced with numbers.")

	    ("huge", po::value<int>()->default_value(5000), "subtrees with at most huge sampled ones would be a out set")

	    ("min-in", po::value<int>()->default_value(100), "subtrees with at least min-in sampled ones would be a out set")
	    ("location_label,l", po::value<vector<string>>()->multitoken(), "location labels, e.g. Germany nonGermany")
	;
/*
	    ("large", po::value<int>()->default_value(1000), "only to show small children aggregated")
	    ("depth", po::value<int>()->default_value(10), "max depth")
	    ("metadata", po::value<string>(), "metadata file, useful if poc is available")
	    ("poc", po::value<string>(), "pangolin of concern. Number of samples with this pangolin is also reported as output.")
*/

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	try {
		po::notify(vm);
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	}

	StateInOut::names = {vm["location_label"].as<vector<string>>()[0], vm["location_label"].as<vector<string>>()[1]};

	string tree_file_name = vm["in"].as<string>();
	INode phylo = load_tree<INode>(tree_file_name) ;

	vector<string> samples_in_vector = load_list_file(vm["in-samples"].as<string>());
	set<string> samples_in(samples_in_vector.begin(), samples_in_vector.end());

	Stat stat(samples_in, 0, vm["huge"].as<int>(), vm["min-in"].as<int>());
	vector<OutSet> samples_out = stat.run(phylo);

	string fn_template = vm["out-samples"].as<string>();
	size_t wc_loc = fn_template.find('?');
	int all_samples_count = 0, all_samples_in_count = 0, all_sample_poor_count=0, all_sample_poor_in_count=0;
	for (size_t i=0; i<samples_out.size(); i++) {
		assert(wc_loc != string::npos);
		string fn = fn_template;
		fn.replace(wc_loc, 1, to_string(i+1));
		ofstream fo(fn);
		int sample_in_count = 0, sample_poor_in_count = 0, sample_poor_count = 0;
		for (size_t j = 0; j < samples_out[i].samples.size(); j++) {
			string s = samples_out[i].samples[j];
			fo << s << endl;
			all_samples_count++;
			if (stat.sample_location_map[s] == StateInOut::IN)
				sample_in_count++;
			if (samples_out[i].is_sample_poor[j]) {
				sample_poor_count++;
				if (stat.sample_location_map[s] == StateInOut::IN)
					sample_poor_in_count++;
			}
		}
		cerr << "Saved " << fn << " with " << samples_out[i].samples.size() << "(in=" << sample_in_count << ") poor=" << sample_poor_count << "(in=" << sample_poor_in_count << ") samples " << endl;
		all_samples_in_count += sample_in_count;
		all_sample_poor_count += sample_poor_count;
		all_sample_poor_in_count += sample_poor_in_count;
	}

	cerr << "Saved " << all_samples_count << " samples in " << samples_out.size() << "(in=" << all_samples_in_count << ") poor=" << all_sample_poor_count << "(in=" << all_sample_poor_in_count << ") files with template " << fn_template << endl;
	
	/*
	cerr << "Selected nodes: ";
	for (auto const& s : stat.selected_nodes) {
		cerr << s << " ";
	}
	cerr << endl;
*/

	return 0;
}
