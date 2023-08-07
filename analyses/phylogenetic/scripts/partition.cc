#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include "tree.h"
#include "state.h"
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/lzma.hpp>
//#include <boost/program_options.hpp>
//namespace po = boost::program_options;
namespace io = boost::iostreams;

typedef int Data;
typedef Node<StateInOut,Data> INode;
array<string, StateInOut::size> StateInOut::names = {"UK", "nonUK"};

template<typename NODE>
struct NodePartitioner : NodePrinterAbstractClass<NODE> {

	string prefix;
	int max_size, min_size;
	int index = 0;
	vector<reference_wrapper<const NODE>> ignored;

	NodePartitioner(string _prefix = "", int _max_size=0, int _min_size=0) : 
		NodePrinterAbstractClass<NODE>(true, true), 
		prefix(_prefix), max_size(_max_size), min_size(_min_size) {}
	
	void partition(const NODE& n) {
		if (n.sample_size >= max_size) {
			for (auto const & c : n.children) {
				partition(c);
			}
		} else if (n.sample_size >= min_size) {
			index++;
			ofstream fo(prefix + to_string(index) + ".tree");
			fo << "tree STATE_0 = ";
			print(fo, n);
			fo << ";" << endl;
		} else {
			//cerr << "IG " << n.sample_size << " " << n.size << " h=" << n.height << endl;
			ignored.push_back(ref(n));
		}
	}

	void print(ostream& os, const NODE& t) {
		if (t.children.size() > 0) {
			os << "(";
			int printed_node_count = 0;
			for (auto &n : t.children) {
				if (printed_node_count > 0)
					os << ",";
				print(os, n);
				printed_node_count ++;
			}
			os << ")";
		}
		this->print_node_info(os, t);
	}
};

int main(int argc, char* argv[]) {
	//po::options_description desc("Allowed options");
	//desc.add_options()
	//    ("help", "Run sankoff")
	//    ("metadata", po::value<string>(), "metadata file")
	//    ("in", po::value<string>(), "input tree")
	//    ("out", po::value<string>(), "output tree")
	//;
	//po::parsed_options parsed_options = po::command_line_parser(argc, argv)
	//	.options(desc)
	//	.run();
	//po::variables_map vm;
	//po::store(parsed_options, vm);

	//string tree = "gisaid-20210324.tree";
	//string tree = vm["in"].as<string>();
	string tree = argv[1];
	//string prefix = "gisaid-20210324-";
	//string prefix = vm["out"].as<string>();
	string prefix = argv[2];
	//string metadata_file = "metadata-gisaid-20210325.tsv.xz";
	//string metadata_file = vm["metadata"].as<string>();
	string metadata_file = argv[3];

	ifstream file(metadata_file, ios_base::in | ios_base::binary);
	io::filtering_streambuf<io::input> in;
	//in.push(gzip_decompressor());
	in.push(io::lzma_decompressor());
	in.push(file);
	std::istream incoming(&in);

	map<string, Metadata> metadata = load_map(incoming);
	cerr << "Metadat loaded " << metadata_file << " " << metadata.size() << endl;

	INode phylo = load_tree<INode>(tree);

	int removed_count;
	phylo.remove_invalid_children(metadata, removed_count);
	cerr << "invalid children removed" << endl;
	NodePartitioner<INode> np(prefix, 10000, 50);
	cerr << "NodePartitioner created" << endl;
	np.partition(phylo);
	cerr << "NodePartitioner partitined" << endl;

	int ignored_count = 0;
	for (auto &i : np.ignored)
		ignored_count += i.get().sample_size;

	cerr << "removed " << removed_count << " and " << np.index << " files saved with pattern " << prefix << "?.tree" << " ignored=" << np.ignored.size() << " and " << ignored_count << " samples" << endl;
	
	return 0;
}
