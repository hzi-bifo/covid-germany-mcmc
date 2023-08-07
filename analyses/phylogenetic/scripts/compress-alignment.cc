#include <unordered_map>
#include <iostream>
#include <fstream>
#include <vector>
#include "fasta.h"
#include <unordered_map>
#include <algorithm>

using namespace std;

struct IdSequence {
	string id, seq;
	IdSequence(string id="", string seq="") : id(id), seq(seq) {}


	bool operator<(const IdSequence& s) const {
		return seq < s.seq;
	}

	bool operator==(const IdSequence& s) const {
		return seq == s.seq;
	}
};

struct FastaLoaderLoader {

	vector<IdSequence> fasta;

	FastaLoaderLoader() {
	}

	void clear() {
		fasta.clear();
	}

	std::string id(std::string raw_id) {
		return raw_id;
	}

	void process(std::string id, std::string seq) {
		fasta.push_back(IdSequence(id, seq));
		if (fasta.size() % 1000 == 0)
			cerr << "loading " << id << " " << fasta.size() << endl;
	}
};

int main(int argc, char* argv[]) {
	FastaLoaderLoader filter;
	FastaLoader<FastaLoaderLoader> fastaLoader(filter);
	fastaLoader.load(cin);

	filter = fastaLoader.filter;

	cerr << "fasta loaded " << filter.fasta.size() << endl;
	int n = filter.fasta.size();
	
	map<int, vector<vector<IdSequence>::iterator>> hash_to_id_seq_index;
	size_t m = 0;
	string column_fix;
	//vector<vector<IdSequence>::iterator> non_duplicates;
	for (vector<IdSequence>::iterator i = filter.fasta.begin(); i != filter.fasta.end(); i++) {
		string id = i->id, seq=i->seq;
		if (m == 0) {
			m = seq.size();
			column_fix = seq;
		} else {
			if (m != seq.size()) {
				cerr << "Incompatible sizes " << m << " " << id << " " << seq.size() << endl;
			}
			for (size_t j=0; j<seq.size(); j++)
				if (column_fix[j] != seq[j])
					column_fix[j] = '*';

			/*
			std::size_t h1 = std::hash<std::string>{}(seq);
			bool duplicate = false;
			if (hash_to_id_seq_index.find(h1) != hash_to_id_seq_index.end()) {
				for (vector<IdSequence>::iterator idx : hash_to_id_seq_index[h1])
					if (seq == idx->seq) {
						duplicate = true;
						break;
					}
			}
			if (!duplicate) {
				non_duplicates.push_back(i);
				hash_to_id_seq_index[h1].push_back(i);
			}
			*/

		}
	}

	int eq_count = 0;
	for (auto & c : column_fix)
		if (c != '*')
			eq_count++;

	sort( filter.fasta.begin(), filter.fasta.end() );
	filter.fasta.erase( unique( filter.fasta.begin(), filter.fasta.end() ), filter.fasta.end() );



	cerr << "column " << column_fix.size() << " eq:" << eq_count << " m:" << m << " size: " << n << " non-dup: " << filter.fasta.size() << " non-dup:" << ""  << endl;
	//for (auto i : non_duplicates) {
	for (auto i = filter.fasta.begin(); i != filter.fasta.end(); i++) {
		string id = i->id, seq = i->seq;
		cout << ">" << id << endl;
		string seq_out = "";
		//for (size_t j=0; j<seq.size(); j++)
		//	if (column_fix[j] == '*')
		//		seq_out += seq[j];
		seq_out = seq;
		cout << seq_out << endl;
	}

	return 0;
}
