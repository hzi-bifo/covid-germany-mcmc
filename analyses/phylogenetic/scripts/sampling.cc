#include "tree.h"
#include "state.h"
#include <boost/program_options.hpp>
#include "rangetree.h"
#include <vector>
#include <numeric>

namespace po = boost::program_options;

array<string, StateInOut::size> StateInOut::names = {"uk", "nonuk"};

struct Data;
typedef Node<StateInOut, Data> INode;

struct Data {
	int dfs_vis, dfs_fin;
	double depth;
	vector<INode*> parents;
	bool back_bone;
	bool sampled;
	Data(int dfs_vis=-1, int dfs_fin=-1, double depth=-1) : dfs_vis(dfs_vis), dfs_fin(dfs_fin), depth(depth), back_bone(false), sampled(false) {
	}
};

ostream& operator<<(ostream& os, const Data& d) {
	os << "(" << d.dfs_vis <<"/" << d.dfs_fin << " " << d.depth << " ";
	for (auto& p : d.parents)
		os << p->label << " ";
	os << d.back_bone << " " << d.sampled << ")";
	return os;
}


//we assume branch_length represets dissimilarity
struct Sampler {
	double range, range_for_state;
	int removed_internal_count, removed_sample_count;
	Sampler(double _range = 0, double range_for_state = 0) : range(_range), range_for_state(range_for_state), removed_internal_count(0), removed_sample_count(0) {
	}

	RangeTree seg;
	vector<INode*> samples;
	vector<Point> node_dfs_values;

	int dfs_vis, dfs_fin;
	void dfs(INode& n, vector<INode*>& parents, double depth) {
		if (n.isLeaf()) {
			samples.push_back(&n);
		}
		n.data = Data(dfs_vis++, -1, depth + n.branch_length);
		for (size_t i=1; parents.size() >= i; i*=2) {
			n.data.parents.push_back(parents[parents.size() - i]);
		}
		//cerr << "D " << n.label << " " << parents.size() << " " << n.data.parents.size() <<  endl;
		parents.push_back(&n);
		for (auto &c: n.children) {
			dfs(c, parents, depth + n.branch_length);
		}
		parents.pop_back();
		n.data.dfs_fin = dfs_fin++;
		node_dfs_values.push_back(Point(n.data.dfs_vis, n.data.dfs_fin));
	}
	
	INode run(INode& n) {
		// we are more restricted for state samples
		assert(range_for_state <= range);
		dfs_vis = dfs_fin = 0;
		node_dfs_values.clear();
		vector<INode*> parents;
		dfs(n, parents, 0);
		cerr << "dfs called" << endl;
		seg = RangeTree(dfs_vis, dfs_fin, node_dfs_values);
		cerr << "seg created" << endl;

		sort(samples.begin(), samples.end(), [](INode* s1, INode* s2) {return s1->data.depth > s2->data.depth;});
		cerr << "samples sorted" << endl;

		for (auto& s: samples) {
			if (!is_covered(*s)) {
				s->data.sampled = true;
				add(*s);
			}
		}

		return sample(n).first;
	}

	void add(INode& n) {
		//cerr << "+ " << n.label << endl;
		INode* nn = &n;
		while (true) { 
			nn->data.back_bone = true;
			if (nn->data.parents.size() > 0) 
				nn = nn->data.parents[0];
			else
				break;
		}
		//cerr << "  seg. " << n.data.dfs_vis << " " << n.data.dfs_fin << "=" << n.data.depth << endl;
		seg.modify(n.data.dfs_vis, n.data.dfs_fin, n.data.depth);
	}

	bool is_covered(INode& n) {
		INode* p = &n;
		//cerr << " ? " << n.label << " " << n.data;
		for (int l = ((int)p->data.parents.size())-1; p->data.parents.size() > 0 && p->data.parents[0]->data.back_bone == false; ) {
			if (l < (int)p->data.parents.size() && p->data.parents[l]->data.back_bone == false)
				p = p->data.parents[l];
			else
				l--;
		}
		//cerr << " first-non-backbone:" << p->label;
		// selected.vis >= p->data.dfs_vis  && selected.fin < p->data.dfs_fin
		if (p->data.parents.size() > 0) {
			Data& data = p->data.parents[0]->data;
			double d = seg.query(data.dfs_vis, 0, dfs_vis, data.dfs_fin),
				min_dis_to_selected_samples = d - data.depth + n.data.depth - data.depth;
			//cerr << " r="<<(d - data.depth + n.data.depth - data.depth) << " bb:" << p->data.parents[0]->label << " d:" << d << " bb-d: " << data.depth << " my-d:" << n.data.depth << endl;
			if ((n.location == StateInOut::IN and min_dis_to_selected_samples <= range_for_state) or (n.location == StateInOut::OUT and min_dis_to_selected_samples <= range)) {
				//if (n.location == StateInOut::IN) {
				//	cerr << "W IN not selected " << n.location << " " << min_dis_to_selected_samples << " " << range_for_state << " " << range << " " << StateInOut::IN << endl;
				//}
				return true;
			} else {
				return false;
			}
		} else {
			//no sample selected yet
			//cerr << "  cov: first" << endl;
			return false;
		}
	}

	pair<INode, bool> sample(const INode& n) {
		//height, size should be recalculated
		INode r(n.label, n.branch_length, n.annotation, list<INode>(), n.location, 1, 1, n.isLeaf() ? 1 : 0);
		for (auto &c: n.children) {
			pair<INode, bool> c_c_ = sample(c);
			if (c_c_.second == false) {
				if (c.isLeaf())
					removed_sample_count++;
				else
					removed_internal_count++;
				continue;
			}
			INode& c_c = c_c_.first;
			if (c_c.children.size() == 1 && !c.isLeaf()) {
				//cerr << "node removed " << c.label << " " << c.branch_length << endl;
				removed_internal_count++;
				for (auto& cc: c_c.children) {
					cc.branch_length += c.branch_length;
					r.height = max(r.height, cc.height + 1);
					r.size += cc.size;
					//r.children.push_back(move(cc));
					r.sample_size += cc.sample_size;
					r.children.push_back(cc);
				}
			} else {
				//normal:
				r.height = max(r.height, c_c.height + 1);
				r.size += c_c.size;
				r.sample_size += c_c.sample_size;
				//r.children.push_back(move(c_c));
				r.children.push_back(c_c);
			}
		}
		return make_pair(r, r.children.size() > 0 || n.data.sampled == true);
	}


	void write_sample_file(string fn) {
		int sample_in = 0, sample_out = 0, all_in=0, all_out=0;
		ofstream fo(fn);
		for (auto const& s: samples) {
			if (s->location == StateInOut::IN)
				all_in++;
			else
				all_out++;
			if (s->data.sampled) {
				if (s->location == StateInOut::IN)
					sample_in++;
				else
					sample_out++;
				fo << s->label << endl;
			}
		}
		cerr << "Samples in=" << sample_in << "/" << all_in << " out=" << sample_out << "/" << all_out << endl;
	}

};

// Algorithm: 
// It removes the newest sample, then every covered case is removed. A sample is covered
//   if there is a sampled case with desired distance (distance=distance over tree edges' branch lengths)
//   to it. For Germany/nonGermany cases this value is different.

int main(int argc, char* argv[]) {

	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help", "Run sankoff")
	    ("short", po::value<double>(), "Branches of length less than or equal to this value are contracted")
	    ("shortstate", po::value<double>(), "Cover length for state (eg. UK)")
	    ("in", po::value<string>(), "input tree")
	    ("out", po::value<string>(), "output tree")
	    ("samples-out", po::value<string>(), "output file to write sample ids")
	    ("location_label,l", po::value<vector<string>>()->multitoken(), "location labels, e.g. Germany nonGermany")
	    ("print-annotation", po::value<bool>()->default_value(true), "print annotation")
	    ("print-internal-node-label", po::value<bool>()->default_value(true), "print internal node labels")
	    ("ilabel", po::value<bool>()->default_value(false), "override internal node labels")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	try {
		po::notify(vm);
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << "\n";
		return 1;
	}

	StateInOut::names = {vm["location_label"].as<vector<string>>()[0], vm["location_label"].as<vector<string>>()[1]};

	string tree_file_name = vm["in"].as<string>(),
		output = vm["out"].as<string>();
	INode phylo = load_tree<INode>(tree_file_name) ;

	Sampler sampler(vm["short"].as<double>(), vm["shortstate"].as<double>());
	INode contracted_phylo = sampler.run(phylo);

	cerr << "removed " << sampler.removed_internal_count << " internal nodes and " << sampler.removed_sample_count << " samples with short branch lenghts" << endl;
	cerr << "  h=" << phylo.height << ">" << contracted_phylo.height << " s=" << phylo.size << ">" << contracted_phylo.size << " ss=" << phylo.sample_size << ">" << contracted_phylo.sample_size << endl;

	if (vm["ilabel"].as<bool>() == true) {
		InternalNodeLabeler<INode> internalNodeLabeler;
		internalNodeLabeler.run(contracted_phylo);
		cerr << "internal nodes relabeled" << endl;
	}
	
	ofstream fo(output);
	NodePrinterGeneral<INode> np(vm["print-annotation"].as<bool>(), vm["print-internal-node-label"].as<bool>());
	//np.print(fo, phylo) << ";" << endl;
	np.print(fo, contracted_phylo) << ";" << endl;

	if (vm.count("samples-out") > 0) {
		sampler.write_sample_file(vm["samples-out"].as<string>());
	}

	cerr << "output saved on " << output<< " " << "[" << vm["print-annotation"].as<bool>() << " " << vm["print-internal-node-label"].as<bool>() << "]" << endl;
	return 0;
}
