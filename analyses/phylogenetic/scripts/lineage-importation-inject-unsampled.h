#include "tree.h"
#include "state.h"
#include <tr1/unordered_map>
#include <iomanip>
#include <string>
#include <set>

using namespace std;


template<typename T>
T to(string v) {
	assert(false);
	return 0;
}

template<>
double to<double>(string v) {
	return stod(v);
}

template<>
string to<string>(string v) {
	return v;
}

string years_to_calendar(double years) {
	tm date_tm = {0};
	date_tm.tm_year = int(years) - 1900;
	date_tm.tm_mday = int( 366 * (years-int(years))) + 1;
	mktime(&date_tm);  // t is now your desired time_t
	char buff[100];
	strftime(buff, sizeof(buff), "%Y-%m-%d", &date_tm);
	return buff;
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

struct Data;
typedef Node<StateInOut, Data> INode;

struct Data {
	double height, year;
	map<string, string> annot_vars;

	int sample_count_in, sample_count_out;

	template<typename T>
	pair<T, T> get_annotation_pair(string key) {
		if ((key == "height_range" || key == "height_95%_HPD") && (annot_vars.find(key) == annot_vars.end() || annot_vars.find(key)->second == "")) {
			assert(annot_vars.find("height") != annot_vars.end());
			T h = to<T>(annot_vars.find("height")->second);
			return make_pair(h, h);
		}
		if (annot_vars.find(key) == annot_vars.end()) {
			cerr << "annot_vars not found " << key << endl;
		}
		assert(annot_vars.find(key) != annot_vars.end());
		string val = annot_vars[key];
		if (!(val[0] == '{' && val[val.size() - 1] == '}')) {
			cerr << "annot_pair  " << val << endl;
		}
		assert(val[0] == '{' && val[val.size() - 1] == '}');
		vector<string> spl = split(val, ',');
		return make_pair(to<T>(spl[0].substr(1)), to<T>(spl[1].substr(0, spl[1].size()-1)));
	}

	double get_location_in_prob() {
        // cerr << "annot_vars:" << annot_vars << " " << endl;
		assert(annot_vars.find("location.set") != annot_vars.end() && annot_vars.find("location.set.prob") != annot_vars.end());
		string lset = annot_vars.find("location.set")->second,
		       lprob = annot_vars.find("location.set.prob")->second;
		if (lset.find(',') != string::npos) {
			pair<double, double> prob_ = get_annotation_pair<double>("location.set.prob");
			// cerr << "prob_ " << prob_.first << " " << prob_.second << endl;
			pair<string, string> lset_ = get_annotation_pair<string>("location.set");
			// cerr << "lset__ " << lset_.first << " " << lset_.second << endl;
			double prob[2] = {prob_.first, prob_.second};
			string lset[2] = {lset_.first, lset_.second};
			for (size_t i=0; i<2; i++) {
                // cerr << "     " << lset[i] << " " <<  StateInOut::names[0] << endl;
				if (lset[i] == StateInOut::names[0])
					return prob[i];
			}
			return 0;
		} else {
			if (lset == "{" + StateInOut::names[0] + "}")
				return 1;
			return 0;
		}
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

template<typename T>
void sort_unique(vector<T>& v) {
	sort(v.begin(), v.end());
	v.erase(unique(v.begin(), v.end()), v.end());
}

//implements a directed distance, loads from mesh output file
struct SampleDistance {
	map<pair<string, string>, double> dist_;
	void load(string fn) {
		ifstream fi(fn);
		string a, b, ratio;
		double pv, x;
		while (fi >> a >> b >> pv >> x >> ratio) {
			dist_[make_pair(a, b)] = pv;
		}
	}

	double dist(string a, string b) {
		if(dist_.find(make_pair(a, b)) == dist_.end()) {
			//cerr << "Dist not found " << a << " " << b << endl;
			return 1;
		}
		assert(dist_.find(make_pair(a, b)) != dist_.end());
		return dist_[make_pair(a, b)];
	}

	vector<string> get_query_samples() {
		vector<string> query_samples;
		for (auto i: dist_) {
			query_samples.push_back(i.first.second);
		}
		//sort(query_samples.begin(), query_samples.end());
		//query_samples.erase(unique(query_samples.begin(), query_samples.end()), query_samples.end());
		sort_unique(query_samples);
		return query_samples;
	}

	/*
	map<string, int> sample_index, query_index;
	vector<string> samples, queries;
	vector<vector<double>> dist_compressed_;
	void compress() {
		for (auto& i : dist_) {
			samples.push_back(i.first.first);
			queries.push_back(i.first.second);
		}
		sort_unique(samples);
		sort_unique(queries);
		for (size_t i = 0; i < samples.size(); i++) {
			sample_index[samples[i]] = i;
		}
		for (size_t i = 0; i < queries.size(); i++) {
			query_index[queries[i]] = i;
		}
		dist_compressed_ = vector<vector<double>> (sample_index.size(), vector<double>(query_index.size(), 1e10));
		for (auto& i : dist_) {
			dist_compressed_[sample_index[i.first.first]][query_index[i.first.second]] = i.second;
		}
	}

	double dist_compressed(size_t a, size_t b) {
		assert(0 <= a && a < samples.size());
		assert(0 <= b && b < queries.size());
		return dist_compressed_[a][b];
	}
	*/

};

struct FullMetadatas {
	vector<string> header;
	map<string, vector<string>> data;
	int header_id_index = -1;

	void load(string fn) {
		ifstream fi(fn);
		for (string line; getline(fi, line); ) {
			vector<string> x = split(line, '\t');
			if (header.size() == 0) {
				header = x;
				for (size_t i=0; i<header.size(); i++) {
					if (header[i] == "Accession.ID")
						header_id_index = i;
				}
			} else {
				assert(header_id_index >= 0);
				data[x[header_id_index]] = x;
			}
		}
	}

	map<string, string> get(string id) const {
		map<string, string> r;
		const vector<string> & val = data.find(id)->second;
		for (size_t i=0; i<header.size(); i++) {
			r[header[i]] = val[i];
		}
		return r;
	}
};

struct Lineage {
	string name;
	vector<string> samples;
	vector<string> query_samples;
	double tmrca = 1e10;
	//string tmrca_annotation = "";
	INode* tmrca_node;
	//If true, the sample for this lineage is not real and should be removed.
	bool allow_new_lineages = false;

	INode* root;
	void set_root(INode& _root) {
		root = &_root;
	}


	//vector<int> samples_compressed;
	//vector<int> query_samples_compressed;

	void add_node(INode& n, const map<string, Metadata>& metadata) {
		if (n.isLeaf()) {
			samples.push_back(n.label);
		}
		if (n.data.year < 2019) {
			cerr << "Invalid year " << n.label << " " << n.data.year << endl;
		}
		//assert(n.data.year >= 2019);
		if (n.data.year < tmrca) {
			tmrca = min(n.data.year, tmrca);
			//tmrca_annotation = n.annotation;
			tmrca_node = &n;
		}
	}

	void set_name(string name) {
		this->name = name;
	}

	double dist(string q, double q_year, SampleDistance& sample_distance) {
		double min_dist = 1e10;
		//cerr << "dist: " << name << " " << q << " " << tmrca << " " << q_year << endl;
		if (q_year >= tmrca) {
			for (auto & s : samples) {
				double sd = sample_distance.dist(s, q);
				if (sd < min_dist) {
					min_dist = sd;
				}
			}
		}
		return min_dist;
	}

	void add_query(string q) {
		query_samples.push_back(q);
	}

	map<string, vector<string>> query_of_samples;
	void add_query(string q, string sample) {
		query_of_samples[sample].push_back(q);
	}

	static void save_cluster_header(ofstream& of) {
		of << "cluster\tseqs\ttmrca\ttmrca_calendar\toldest\tmostrecent\ttreefile\tcutoff\ttmrca_HPD_upper\ttmrca_HPD_lower\ttmrca_range_upper\ttmrca_range_lower"<< endl;
	}

	void save_cluster(ofstream& of, const map<string, Metadata>& metadata, double ref, string treefile, ofstream& of_singles, string cutoff) {
		double most_recent = -1e10, oldest = 1e10;
		for (string s : samples) {
			double d = years(metadata.find(s)->second.date);
			most_recent = max(most_recent, d);
			oldest = min(oldest, d);
		}
		for (string s : query_samples) {
			double d = years(metadata.find(s)->second.date);
			most_recent = max(most_recent, d);
			oldest = min(oldest, d);
		}
		//add, check why it was not here?
		int query_of_samples_size = 0;
		for (auto sq : query_of_samples) {
			query_of_samples_size += sq.second.size();
			string s = sq.first;
			double d = years(metadata.find(s)->second.date);
			most_recent = max(most_recent, d);
			oldest = min(oldest, d);
		}
		//map<string, string> annot_vars = get_tmrca_annotation_var();
		//pair<double, double> tmrca_HPD = annot_pair(annot_vars["height_95%_HPD"], annot_vars),
		//		tmrca_range =  annot_pair(annot_vars["height_range"], annot_vars);
		//map<string, string> annot_vars = get_tmrca_annotation_var();
		//cerr << "Saving cluster " << name << " " << tmrca_node->annotation << endl;
		pair<double, double> tmrca_HPD = tmrca_node->data.get_annotation_pair<double>("height_95%_HPD"),
				tmrca_range =  tmrca_node->data.get_annotation_pair<double>("height_range");
		if (samples.size() + query_samples.size() + query_of_samples_size == 1) {
			of_singles << name << "\t" << (samples.size() + query_samples.size() + query_of_samples_size) << "\t" << tmrca << "\t" << years_to_calendar(tmrca) << "\t" << oldest << "\t" << most_recent << "\t" << treefile << "\t" 
				<< cutoff << "\t" << tmrca_HPD.second + ref << "\t" << tmrca_HPD.first +ref << "\t" << tmrca_range.second +ref << "\t" << tmrca_range.first +ref << endl;
		} else {
			of << name << "\t" << (samples.size() + query_samples.size() + query_of_samples_size) << "\t" << tmrca << "\t" << years_to_calendar(tmrca) << "\t" << oldest << "\t" << most_recent << "\t" << treefile << "\t" 
				<< cutoff << "\t" << tmrca_HPD.second + ref << "\t" << tmrca_HPD.first +ref << "\t" << tmrca_range.second +ref << "\t" << tmrca_range.first +ref << endl;
		}
	}

	static vector<string> cluster_samples_header_names;
	static void save_cluster_samples_header(ofstream& fo) {
		bool first = true;
		for (auto& h : cluster_samples_header_names) {
			if (!first) fo << "\t";
			fo << h;
			first = false;
		}
		fo << endl;
	}

	void save_cluster_sample(ofstream& fo, string s, const FullMetadatas& full_metadatas) {
		//cerr << "save_cluster_sample " << s << endl;
		map<string, string> mv = full_metadatas.get(s);
		mv["sample_date"] = mv["Collection.date"];
		mv["decimal_date"] = to_string(years(mv["sample_date"]));
		mv["taxon_label"] = mv["Accession.ID"];
		mv["is_sampled"] = "FALSE";
		mv["cluster"] = name;


		bool first = true;
		for (auto& h : cluster_samples_header_names) {
			if (!first) fo << "\t";
			fo << mv[h];
			first = false;
		}
		fo << endl;
	}

	void save_cluster_samples(ofstream& of, const FullMetadatas& full_metadatas) {
		//cerr << "save_cluster_samples " << samples << ", " << query_samples << endl;
		for (string s : samples) {
			save_cluster_sample(of, s, full_metadatas);
		}
		for (string s : query_samples) {
			save_cluster_sample(of, s, full_metadatas);
		}
		// added here, why was not here?
		for (auto qs : query_of_samples) {
			for (auto q : qs.second)
				save_cluster_sample(of, q, full_metadatas);
		}
	}

	/*
	void compress(SampleDistance& sample_distance) {
		for (auto s : samples)
			samples_compressed.push_back(sample_distance.sample_index[s]);
		for (auto s : query_samples)
			query_samples_compressed.push_back(sample_distance.sample_index[s]);
	}

	double dist_compressed(int q, double q_year, SampleDistance& sample_distance) {
		double min_dist = 1e10;
		//cerr << "dist: " << name << " " << q << " " << tmrca << " " << q_year << endl;
		if (q_year >= tmrca) {
			for (auto & s : samples_compressed) {
				double sd = sample_distance.dist_compressed(s, q);
				if (sd < min_dist) {
					min_dist = sd;
				}
			}
		}
		return min_dist;
	}
	*/
};
vector<string> Lineage::cluster_samples_header_names = {"Virus.name", "Type", "Accession.ID", "Collection.date", "Location", "Additional.location.information", "Sequence.length", "Host", "Patient.age", "Gender", "Clade", "Pango.lineage", "Pangolin.version", "Variant", "AA.Substitutions", "Submission.date", "Is.reference.", "Is.complete.", "Is.high.coverage.", "Is.low.coverage.", "N.Content", "GC.Content", "sample_date", "decimal_date", "taxon_label", "is_sampled", "cluster"};

struct TreeLeafLabelAssignerIdentical {
    string operator()(string original_label) {
        return original_label;
    }
};

struct Lineages {
	list<Lineage> lineages_NA, lineages_half;

	const map<string, Metadata>& metadata;

	Lineages(const map<string, Metadata>& metadata) : metadata(metadata) {
	}

	void save_clusters(string fn, double reference_height, string treefile, string fn_singles, list<Lineage>& lineages, string cutoff = "NA") {
		ofstream fo(fn);
		ofstream fo_singles(fn_singles);
		fo << std::fixed << setprecision(4);
		fo_singles << std::fixed << setprecision(4);

		Lineage::save_cluster_header(fo);
		Lineage::save_cluster_header(fo_singles);
		for (auto& l: lineages) {
			l.save_cluster(fo, metadata, reference_height, treefile, fo_singles, cutoff);
		}
	}

	void save_cluster_samples(string fn, const FullMetadatas& full_metadatas, list<Lineage>& lineages) {
		ofstream fo(fn);
		fo << std::fixed << setprecision(4);
		Lineage::save_cluster_samples_header(fo);
		for (auto& l: lineages) {
			l.save_cluster_samples(fo, full_metadatas);
		}
	}

	template<typename TreeLeafLabelAssigner>
	pair<string, double> lineage_tree(const FullMetadatas& full_metadatas, Lineage& l, INode& n, double loc_p_cutoff, TreeLeafLabelAssigner& tree_leaf_label_assigner) {
		StateInOut loc = n.location;
		if (loc_p_cutoff != -1) {
			if (n.data.get_location_in_prob() >= loc_p_cutoff)
				loc = StateInOut::IN;
			else
				loc = StateInOut::OUT;
		}
		if (loc == StateInOut::OUT) 
			return make_pair("", 0.0);
		string ret = "";
		int child_count = 0;
		double prev_branch_length = 0;
		for (INode& c : n.children) {
			pair<string, double> r = lineage_tree<TreeLeafLabelAssigner>(full_metadatas, l, c, loc_p_cutoff, tree_leaf_label_assigner);
			if (r.first != "") {
				// the node have more than one child, print previous length and act normally
				if (child_count == 1)
					ret += ":" + std::to_string(prev_branch_length);
				if (child_count > 0) {
					ret += ",";
				}
				ret += r.first;
				if (child_count > 0) {
					ret += ":" + std::to_string(r.second);
				}
				prev_branch_length = r.second;
				child_count++;
			}
		}
		if (child_count > 1) {
			ret = "(" + ret + ")";
		}
		if (n.isLeaf()) {
			ret += tree_leaf_label_assigner(n.label);
			if (l.query_of_samples.find(n.label) != l.query_of_samples.end()) {
				ret += ":" + std::to_string(0.0);
				for (auto & i : l.query_of_samples.find(n.label)->second) {
					ret += "," + i + ":" + std::to_string(0.0);
				}
				ret = "(" + ret + ")";
			}
		}
		//ret += ":" + to_string(n.branch_length);
		return make_pair(ret, n.branch_length + (child_count > 1 ? 0.0 : prev_branch_length));
	}

	template<typename TreeLeafLabelAssigner>
	void save_lineage_trees(ostream& fo, const FullMetadatas& full_metadatas, list<Lineage>& lineages, double loc_p_cutoff, TreeLeafLabelAssigner& tree_leaf_label_assigner) {
		//ofstream fo(fn);
		//fo << std::fixed << setprecision(4);
		for (auto& l: lineages) {
			fo << l.name << "=" << lineage_tree(full_metadatas, l, *l.root, loc_p_cutoff, tree_leaf_label_assigner).first << endl;
		}
	}

};

struct LineageExtraction {

	const map<string, Metadata>& metadata; 

	double reference_height = -1;

	Lineages& lin;

	LineageExtraction(const map<string, Metadata>& metadata, Lineages& lin) :
		metadata(metadata), lin(lin) {
	}


	void dfs(INode& n1, double h = 0) {
		n1.data.height = h;
		n1.data.set_annotation_vars(n1.annotation);
		if (n1.isLeaf()) {
			if (!(metadata.find(n1.label) != metadata.end())) {
				cerr << "L not in metadata " << n1.label << endl;
			}
			assert(metadata.find(n1.label) != metadata.end());
			if (reference_height == -1) {
				reference_height = years(metadata.find(n1.label)->second.date) - h;
			} else {
				if (!(abs(reference_height - (years(metadata.find(n1.label)->second.date) - h)) < 0.01)) {
					cerr << "Incompatible heights: " << n1.label << " metadata=" << years(metadata.find(n1.label)->second.date) << " h=" << h << " " << reference_height << endl;
				}
				assert(abs(reference_height - (years(metadata.find(n1.label)->second.date) - h)) < 0.01);
			}
		}
		//cerr << "D " << n1.label << " " << n2.label << endl;
		for (auto i1 = n1.children.begin(); i1 != n1.children.end(); i1++) {
			dfs(*i1, h + i1->branch_length);
		}
		assert (reference_height != -1);
		n1.data.year = reference_height + h;
	}

	void assign_date(INode& n) {
		dfs(n, 0);
	}

	void find_lineages(INode& n, Lineage* current_lineage, string name_prefix, list<Lineage>& lineages, double loc_p_cutoff, bool allow_new_lineages, ostream& log_file) {
		StateInOut loc = n.location;
		if (loc_p_cutoff != -1) {
			//cerr << "AP " << n.annotation << " " << n.data.annot_vars["location.set.prob"] << endl;
			/*
			pair<double, double> prob_ = n.data.get_annotation_pair<double>("location.set.prob");
			cerr << "prob_ " << prob_.first << " " << prob_.second << endl;
			pair<string, string> lset_ = n.data.get_annotation_pair<string>("location.set");
			cerr << "lset__ " << lset_.first << " " << lset_.second << endl;
			double prob[2] = {prob_.first, prob_.second};
			string lset[2] = {lset_.first, lset_.second};
			// if location is not in location.set, prob=0
			if (loc_p_cutoff == 0)
				loc = StateInOut::IN;
			else
				loc = StateInOut::OUT;
			for (size_t i=0; i<2; i++) {
				if (lset[i] == StateInOut::names[0])
					if (prob[i] > loc_p_cutoff)
						loc = StateInOut::IN;
			}
			*/
			if (n.data.get_location_in_prob() >= loc_p_cutoff)
				loc = StateInOut::IN;
			else
				loc = StateInOut::OUT;
		}
        // cerr << "find_lineages " << n.label << " " << loc << " " << n.data.get_location_in_prob() << endl;
		//if (n.location != StateInOut::IN) {
		if (loc != StateInOut::IN) {
			//if (n.isLeaf()) {
			//	cerr << "sample out " << n.location << endl;
			//}
			current_lineage = 0;
		} else {
			if (current_lineage == 0) {
				lineages.push_back(Lineage());
				current_lineage = &lineages.back();
				current_lineage->set_name(name_prefix + to_string(lineages.size()));
				current_lineage->set_root(n);
				if (log_file) {
					log_file << "New lineage" << name_prefix + to_string(lineages.size()) << " sub-tree sample count (in/out) " << n.data.sample_count_in << " " << n.data.sample_count_out << " " << endl;
				}
			}
			current_lineage->add_node(n, metadata);
		}
		// We create a lineage for each of the StateInOut::OUT samples,
		//   then when queries are assigned to lienages, in this case (allow_new_lineages), new lineages are created!
		if (allow_new_lineages && loc != StateInOut::IN && n.isLeaf()) {
			lineages.push_back(Lineage());
			lineages.back().add_node(n, metadata);
			lineages.back().allow_new_lineages = true;
			lineages.back().set_name(name_prefix + to_string(lineages.size()));
		}
		for (INode& c : n.children) {
			find_lineages(c, current_lineage, name_prefix, lineages, loc_p_cutoff, allow_new_lineages, log_file);
		}
		/*
		if (current_lineage == &lineage) {
			lineages.push_back(lineage);
		}*/
	}
};

struct QueryAssignment {

	Lineages& lin;

	QueryAssignment(Lineages& lin) : lin(lin) {
	}

	/*
	//assign query samples (unsampled) to lineages
	void assign_to_lineages(const vector<string>& query_samples, SampleDistance& sample_distance, bool compressing = false) {
		if (compressing) {
			sample_distance.compress();
			for (auto &lin: lineages) {
				lin.compress(sample_distance);
			}
		}
		vector<string> not_assigned_queries;
		for (const string& q: query_samples) {
			double q_year = years(metadata.find(q)->second.date);
			double min_dist = 1e10;
			Lineage* min_lin = 0;
			if (compressing) {
				int q_index = sample_distance.query_index[q];
				for (auto &lin: lineages) {
					if (lin.dist_compressed(q_index, q_year, sample_distance) < min_dist)
						min_lin = &lin;
				}
			} else {
				for (auto &lin: lineages) {
					if (lin.dist(q, q_year, sample_distance) < min_dist)
						min_lin = &lin;
				}
			}
			if (min_lin != 0) {
				min_lin->add_query(q);
			} else {
				not_assigned_queries.push_back(q);
			}
		}
		if (not_assigned_queries.size() > 0) {
			cerr << "W not lineage for queries (" << not_assigned_queries.size() << ") ";
			int cnt = 0;
			for (string q : not_assigned_queries) {
				cerr << q << " ";
				cnt++;
				if (cnt == 2) break;
			}
			cerr << endl;
		}
	}
	*/


	//when distance file is so large!
	void assign_to_lineages(string sample_distance_fn, bool allow_new_lineages, double dist_threshold, const map<string, Metadata>& metadata, ostream& log_file) {
		map<string, Lineage*> sample_lineage_NA_map, sample_lineage_half_map;
		for (Lineage& l : lin.lineages_NA) {
			for (string s: l.samples) {
				sample_lineage_NA_map[s] = &l;
			}
		}
		for (Lineage& l : lin.lineages_half) {
			for (string s: l.samples) {
				sample_lineage_half_map[s] = &l;
			}
		}
		cerr << "sample_lineage_map assigned " << sample_lineage_NA_map.size() << " " << sample_lineage_half_map.size() << endl;
		set<string> all_queries; // to track queries which did not assign to any lineage
		map<string, tuple<double, Lineage*, string>> query_min_dist_lineage_NA, query_min_dist_lineage_half;
		string s, q, ratio;
		double pv, dist;
		ifstream fi(sample_distance_fn);
		cerr << "loading " << sample_distance_fn << endl;
		int rejection_dist_threshold_cnt = 0, rejection_dist_threshold_not_cnt = 0;
		for (int cnt=0; fi >> s >> q >> dist >> pv >> ratio; cnt++) {
			all_queries.insert(q);
			if (cnt % 10000 == 0) {
				cerr << "\r sample distance line: " << cnt << " " << s << " " << q << " " << pv << " " << dist << " ratio=" << ratio << " " << " qmd=" << query_min_dist_lineage_NA.size() << " " << query_min_dist_lineage_half.size() << " dtr:" << rejection_dist_threshold_cnt << "+" << rejection_dist_threshold_not_cnt << "      " << flush;
			}
			if (dist > dist_threshold) {
				rejection_dist_threshold_cnt++;
				continue;
			} else {
				rejection_dist_threshold_not_cnt++;
			}
			//cerr << "reading " << cnt << " " << s << endl;
			//dist_[make_pair(a, b)] = pv;
			auto q_i_NA = query_min_dist_lineage_NA.find(q), q_i_half = query_min_dist_lineage_half.find(q);
			//double md_NA = query_min_dist_lineage.find(q) != query_min_dist_lineage.end() ? query_min_dist_lineage.find(q)->second.first : 1e10;
			double md_NA = q_i_NA != query_min_dist_lineage_NA.end() ? get<0>(q_i_NA->second) : 1e10,
			       md_half = q_i_half != query_min_dist_lineage_half.end() ? get<0>(q_i_half->second) : 1e10;
			double q_year = years(metadata.find(q)->second.date);


			auto slm_NA_i = sample_lineage_NA_map.find(s), 
			     slm_half_i = sample_lineage_half_map.find(s);
			//DEBUG
			//if (sample_lineage_map[s]->tmrca < 2019 or q_year < 2019) {
			//	cerr << "Date problem " << std::fixed << setprecision(4) << s << " " << q << " tmrca=" << sample_lineage_map[s]->tmrca << " q-year=" << q_year << " pv=" << pv << " md=" << md << endl;
			//}
			//END DEBUG
			if (slm_NA_i != sample_lineage_NA_map.end()) {
				if (slm_NA_i->second->tmrca <= q_year + 1.0/365 && dist < md_NA) {
					query_min_dist_lineage_NA[q] = make_tuple(dist, slm_NA_i->second, s);
				}
			}
			if (slm_half_i != sample_lineage_half_map.end()) {
				if (slm_half_i->second->tmrca <= q_year + 1.0/365 && dist < md_half) {
					query_min_dist_lineage_half[q] = make_tuple(dist, slm_half_i->second, s);
				}
			}
			/*
			if (sample_lineage_map.find(s) == sample_lineage_map.end())
				continue;
			assert(sample_lineage_map[s]->tmrca > 2019 and q_year > 2019);
			if (sample_lineage_map[s]->tmrca <= q_year + 1.0/365 && dist < md) {
				query_min_dist_lineage[q] = make_pair(dist, sample_lineage_map[s]);
			} else {
				//cerr << "TMRCA problem " << std::fixed << setprecision(4) << s << " " << q << " tmrca=" << sample_lineage_map[s]->tmrca << " q-year=" << q_year << " pv=" << pv << " md=" << md << endl;
			}*/
		}
		cerr << endl;
		cerr << "sample distance loaded" << query_min_dist_lineage_NA.size() << " " << query_min_dist_lineage_half.size() << " among " << all_queries.size() << " " << "(not assigned " << (all_queries.size() - query_min_dist_lineage_NA.size()) << " " << (all_queries.size() - query_min_dist_lineage_half.size()) << ", threshold=" << dist_threshold << ")" << endl;
		//find non-injected queries ???
		if (log_file)
			log_file << "query_min_dist_lineage_NA (lineage query sample distance)" << endl;
		for (auto& i : query_min_dist_lineage_NA) {
			get<1>(i.second)->add_query(i.first, get<2>(i.second));
			if (log_file) {
				log_file << "A " << get<1>(i.second)->name << " " << i.first << " " << get<2>(i.second) << " " << get<0>(i.second) << " " << endl;
			}
		}
		if (log_file)
			log_file << "query_min_dist_lineage_half (lineage query sample distance)" << endl;
		for (auto& i : query_min_dist_lineage_half) {
			get<1>(i.second)->add_query(i.first, get<2>(i.second));
			if (log_file) {
				log_file << "A " << get<1>(i.second)->name << " " << i.first << " " << get<2>(i.second) << " " << get<0>(i.second) << " " << endl;
			}
		}
		cerr << "lin assigned " << query_min_dist_lineage_NA.size() << " " << query_min_dist_lineage_half.size() << endl;
		if (log_file) {
			log_file << "lin assigned " << query_min_dist_lineage_NA.size() << " " << query_min_dist_lineage_half.size() << endl;
		}
	}

	// If allow_new_lineages == true, we add lineages with samples.location==OUT,
	//  this function removes these lineages
	void remove_empty_lineages() {
		array<list<Lineage>*, 2> lineages_containers{&lin.lineages_NA, &lin.lineages_half};
		for (auto & lineages : lineages_containers) {
			for (auto i = lineages->begin(); i != lineages->end(); i++) {
				//These lineages should have exactly one sample from OUT
				if (i->allow_new_lineages) {
					assert(i->samples.size() == 1);
					i->samples.clear();
				}
			}
			auto end = remove_if(lineages->begin(), lineages->end(),
				  [](const Lineage& l) {
				  	return l.samples.size() + l.query_samples.size() == 0;
				});
			lineages->erase(end, lineages->end());
		}
	}


};

//when we only want to cluster query sequences based on their distance, the tree is ignored
struct QueryClustering {

	Lineages& lin;
	string name_prefix;

	set<string> filter_out;

	QueryClustering(Lineages& lin, string name_prefix) : lin(lin), name_prefix(name_prefix) {
	}

	void load_filter_out(string fn) {
		ifstream fi(fn);
		for (string l; fi >> l; )
			filter_out.insert(l);
		cerr << "filter-out loaded " << fn << " " << filter_out.size() << endl;
	}

	string min_string(string a, string b, string c) {
		if (a != "")
			return min(a, min(b, c));
		return min(b, c);
	}

	string find_root(string q, map<string, string>& query_min_cluster_mate) {
		auto i = query_min_cluster_mate.find(q);
		if (i == query_min_cluster_mate.end() || q == i->second)
			return q;
		return query_min_cluster_mate[q] = find_root(i->second, query_min_cluster_mate);
	}

	void assign_to_lineages(string sample_distance_fn, double dist_threshold, const map<string, Metadata>& metadata) {
		set<string> all_queries; // to track queries which did not assign to any lineage
		map<string, string> query_min_cluster_mate;
		string s, q, ratio;
		double pv, dist;
		ifstream fi(sample_distance_fn);
		cerr << "loading " << sample_distance_fn << endl;
		int rejection_dist_threshold_cnt = 0, rejection_dist_threshold_not_cnt = 0;
		for (int cnt=0; fi >> s >> q >> dist >> pv >> ratio; cnt++) {
			if (filter_out.find(q) != filter_out.end() || filter_out.find(s) != filter_out.end())
				continue;
			all_queries.insert(q);
			all_queries.insert(s);
			if (cnt % 10000 == 0) {
				cerr << "\r sample distance line: " << cnt << " " << s << " " << q << " " << pv << " " << dist << " ratio=" << ratio << " " << " qmd=" << query_min_cluster_mate.size() << " dtr:" << rejection_dist_threshold_cnt << "+" << rejection_dist_threshold_not_cnt << "      " << flush;
			}
			if (dist > dist_threshold) {
				rejection_dist_threshold_cnt++;
				continue;
			} else {
				rejection_dist_threshold_not_cnt++;
			}

			query_min_cluster_mate[q] = min_string(query_min_cluster_mate.find(q) != query_min_cluster_mate.end() ? query_min_cluster_mate[q] : "", q, s);
			query_min_cluster_mate[s] = min_string(query_min_cluster_mate.find(s) != query_min_cluster_mate.end() ? query_min_cluster_mate[s] : "", q, s);
		}
		cerr << endl;

		cerr << "sample distance loaded" << " among " << all_queries.size() << " " << "(not assigned " << (all_queries.size() - query_min_cluster_mate.size()) << " " << ", threshold=" << dist_threshold << ")" << endl;

		map<string, Lineage*> root_lineage;
		for (auto & qmcm : query_min_cluster_mate) {
			string r = find_root(qmcm.first, query_min_cluster_mate);
			if (root_lineage.find(r) == root_lineage.end()) {
				lin.lineages_NA.push_back(Lineage());
				lin.lineages_NA.back().set_name(name_prefix + to_string(lin.lineages_NA.size()));
				root_lineage[r] = &lin.lineages_NA.back();
			}
			root_lineage[r]->add_query(qmcm.first);
		}

		cerr << "lin assigned " << query_min_cluster_mate.size() << endl;
	}
};

void fill_tree_info(INode& n, const FullMetadatas& metadata) {
	if (n.isLeaf()) {
		bool is_in = n.location == StateInOut::IN;
		n.data.sample_count_in = is_in;
		n.data.sample_count_out = !is_in;
	} else {
		n.data.sample_count_in  = 0;
		n.data.sample_count_out = 0;
		for (INode& c : n.children) {
			fill_tree_info(c, metadata);
			n.data.sample_count_in  += c.data.sample_count_in;
			n.data.sample_count_out += c.data.sample_count_out;
		}
	}
}
