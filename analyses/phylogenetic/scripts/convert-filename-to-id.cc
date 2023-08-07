#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <cassert>
using namespace std;

string trim(string s) {
	s.erase(0, s.find_first_not_of("\t\n\v\f\r ")); // left trim
	s.erase(s.find_last_not_of("\t\n\v\f\r ") + 1); // right trim
	return s;
}



vector<string> split(string to_split, string delimiter) {
	size_t pos = 0;
	vector<string> matches{};
	do {
		pos = to_split.find(delimiter);
		int change_end;
		if (pos == string::npos) {
			pos = to_split.length() - 1;
			change_end = 1;
		}
		else {
			change_end = 0;
		}
		matches.push_back(to_split.substr(0, pos+change_end));

		to_split.erase(0, pos+1);

	}
	while (!to_split.empty());
	return matches;

}

map<string, string> filename_id;
string get_filename_id(string fn) {
	auto i = filename_id.find(fn);
	if (i == filename_id.end()) {
		ifstream fi(fn);
		string line; getline(fi, line);
		//line = trim(line);
		assert(line.size() > 0 && line[0] == '>');
		line = line.substr(1);
		if (line.find('|') != string::npos) {
			line = line.substr(0, line.find('|'));
		}
		filename_id[fn] = line;
		i = filename_id.find(fn);
	}
	//return filename_id[fn];
	return i->second;
}

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
	int val = (int) (percentage * 100);
	int lpad = (int) (percentage * PBWIDTH);
	int rpad = PBWIDTH - lpad;
	printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
	fflush(stdout);
}

int main() {
	long long cnt = 0;
	for (char line[1000]; cin.getline(line, sizeof(line)); cnt++) {
		int index = 0;
		for (char* c = line, *last_tab=line; *c; c++) {
			if (*c == '\t' && index < 2) {
				*c = 0;
				cout << get_filename_id(last_tab) << "\t";
				index++;
				last_tab = c+1;
				if (index == 2) {
					cout << last_tab << endl;
					break;
				}
			} 
		}
		/*
		auto x = split(line, "\t");
		for (int i=0; i<2; i++)
			x[i] = get_filename_id(x[i]);
		for (size_t i=0; i<x.size(); i++) {
			if (i+1 != x.size())
				cout << "\t";
			cout << x[i];
		}
		cout << endl;
		*/
		if (cnt % 10000 == 0) {
			fprintf(stderr, "\r%15lld", cnt);
			fflush(stderr);
		}
	}
	return 0;
}
