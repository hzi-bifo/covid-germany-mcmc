#include <iomanip>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <assert.h>
using namespace std;

set<string> seqs;

void load_distances(const char* fn) {
		ifstream fi(fn);
		cerr << "loading " << fn << endl;
		string s, q, ratio;
		double pv, dist;
		for (int cnt=0; fi >> s >> q >> dist >> pv >> ratio; cnt++) {
			if (seqs.find(q) != seqs.end() && seqs.find(s) != seqs.end()) {
				cout << s << "\t" << q << "\t" << dist << "\t" << pv << "\t" << ratio << endl;
			}
		}
}

void load_sequences(const char* fn) {
	ifstream fi(fn);
	for (string line; fi >> line; ) {
		seqs.insert(line);
	}
}

int main(int argc, char* argv[]) {
	assert(argc == 3);
	load_sequences(argv[1]);
	load_distances(argv[2]);
	return 0;
}
