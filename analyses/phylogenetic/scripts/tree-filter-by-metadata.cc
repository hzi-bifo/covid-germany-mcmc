#include "tree.h"
#include "state.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

typedef Node<StateInOut, int> INode;

struct SampleInMetadataIncluder {
	const map<string, Metadata>& id_to_name;
	SampleInMetadataIncluder(const map<string, Metadata>& id_to_name) : id_to_name(id_to_name) {}
	bool operator()(const INode& n) {
		return n.isLeaf() ? (id_to_name.find(n.label) != id_to_name.end()) : false;
	}
};

array<string, StateInOut::size> StateInOut::names = {"uk", "nonuk"};

int main(int argc, char* argv[]) {
	po::options_description desc("Allowed options");
	desc.add_options()
	    ("in", po::value<string>(), "input tree")
	    ("out", po::value<string>(), "output tree")
	    ("metadata", po::value<string>(), "metadata file")
	    ("samples", po::value<string>(), "sample files")
	;
	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
		.options(desc)
		.run();
	po::variables_map vm;
	po::store(parsed_options, vm);

	string annotation = vm["metadata"].as<string>();
	map<string, Metadata> id_to_name = load_map(annotation);
	cerr << "metadata loaded" << endl;

	string tree_file_name = vm["in"].as<string>(), output = vm["out"].as<string>();
	INode phylo0 = load_tree<INode>(tree_file_name) ;

	SampleInMetadataIncluder sampleInMetadataIncluder(id_to_name);
	SubtreeExtractorOverSamples<INode, SampleInMetadataIncluder> cleaner(sampleInMetadataIncluder);
	INode phylo2 = cleaner.run(phylo0);

	if (vm.count("samples") > 0) {
		SamplePrinter<INode> sp;
		sp.run(phylo2, vm["samples"].as<string>());
		cerr << "samples printed in " << vm["samples"].as<string>() << " for " << sp.printed_count << " samples" << endl;
	}
	NodePrinterGeneral<INode> np(false, false, false);
	ofstream fo(output);
	np.print(fo, phylo2) << ";" << endl;

	cerr << "done " << tree_file_name << " -> " << output<< " " << "" << endl;
	return 0;
}
