#include "BubbleChain.hpp"
#include "GFAReader.hpp"
#include "boost/program_options.hpp"
#include <unordered_set>
#include <stdexcept>
#include <iostream>

using std::unordered_set;
using std::cout;
using std::ofstream;
using std::getline;
using std::runtime_error;
using std::experimental::filesystem::create_directories;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void parse_nodes_from_assembly_summary_line(string& line, unordered_set<string>& forward_nodes, unordered_set<string>& reverse_nodes){
    ///
    /// Parse a line of the AssemblySummary.csv with the following format:
    ///     Rank,EdgeId,EdgeIdRc,Length,CumulativeLength,LengthFraction,CumulativeFraction
    ///     0,800526,800527,169589,169589,5.80388e-05,5.80388e-05
    ///

    string token;
    uint64_t n_commas = 0;

    for (auto& c: line){
        if (c == ','){
            if (n_commas == 1){
                forward_nodes.insert(token);
            }
            else if (n_commas == 2){
                reverse_nodes.insert(token);
                break;
            }

            token.resize(0);
            n_commas++;
        }
        else{
            token += c;
        }
    }
}


void extract_node_sets_from_assembly_summary(
        path assembly_summary_path,
        unordered_set<string>& forward_nodes,
        unordered_set<string>& reverse_nodes) {

    ifstream assembly_summary(assembly_summary_path);

    if (not assembly_summary.is_open()){
        throw runtime_error("ERROR: could not open bubble chain file: " + assembly_summary_path.string());
    }

    string line;
    uint64_t l = 0;
    while (getline(assembly_summary, line)) {
        if (l == 0) {
            l++;
            continue;
        }

        if (line.empty()) {
            throw runtime_error("ERROR: empty line found in file: " + assembly_summary_path.string() + " at line index " + std::to_string(l));
        }

        parse_nodes_from_assembly_summary_line(line, forward_nodes, reverse_nodes);

        l++;
    }
}


void extract_bubble_chains_from_gfa(path gfa_path, path bubble_path, path assembly_summary_path, path output_dir){
    create_directories(output_dir);
    ifstream bubble_chain_file(bubble_path);

    if (not bubble_chain_file.is_open()){
        throw runtime_error("ERROR: could not open bubble chain file: " + bubble_path.string());
    }

    path output_path = gfa_path.filename();
    output_path.replace_extension("bubble_chains.gfa");
    output_path = output_dir / output_path;

    ofstream output_gfa(output_path);

    if (not output_gfa.good()){
        throw runtime_error("ERROR: output GFA could not be written: " + output_path.string());
    }

    unordered_set<string> forward_nodes;
    unordered_set<string> reverse_nodes;

    extract_node_sets_from_assembly_summary(assembly_summary_path, forward_nodes, reverse_nodes);

    GFAReader gfa_reader(gfa_path);
//    gfa_reader.map_sequences_by_node();

    string line;
    uint64_t l = 0;

    BubbleChainComponent chain_component;
    BubbleChainComponent previous_chain_component;
    uint64_t n_forward = 0;
    uint64_t n_reverse = 0;

    ///
    /// Parse the bubble chain CSV with the format:
    ///     Chain,Circular,Position,Segment0,Segment1,Segment2,Segment3,Segment4,
    ///
    string s;
    while (getline(bubble_chain_file, line)){
        if (l == 0){
            l++;
            continue;
        }
        if (line.empty()){
            throw runtime_error("ERROR: empty line found in file: " + bubble_path.string() + " at line index " + std::to_string(l));
        }

        parse_line_as_bubble_chain_component(line, chain_component);

        if (chain_component.chain != previous_chain_component.chain){
            cout << n_forward << " " << n_reverse << '\n';
            n_forward = 0;
            n_reverse = 0;
        }

        for (auto& segment: chain_component.segments) {
            if (forward_nodes.count(segment) > 0){
                n_forward++;
            }
            else if (reverse_nodes.count(segment) > 0){
                n_reverse++;
            }
            else {
                throw runtime_error("ERROR: node not in set of known forward/reverse nodes from AssemblySummary.csv: " + segment);
            }
        }

//        for (auto& segment: chain_component.segments){
//            gfa_reader.read_line(s, gfa_reader.sequence_line_indexes_by_node.at(segment));
//            output_gfa << s;
//            cout << s << '\n';
//        }


        previous_chain_component = chain_component;
        l++;
    }

}


int main(int argc, char* argv[]){
    path gfa_path;
    path bubble_path;
    path assembly_summary_path;
    path output_dir;

    options_description options("Arguments");

    options.add_options()
            ("gfa",
             value<path>(&gfa_path),
             "File path of GFA file containing shasta assembly graph")

            ("bubbles",
             value<path>(&bubble_path),
             "File path of shasta BubbleChains.csv output which describes the bubble chains contained in the GFA")

            ("summary",
             value<path>(&assembly_summary_path),
             "File path of shasta AssemblySummary.csv output which describes the nodes contained in the GFA")

            ("output_dir",
             value<path>(&output_dir)->
             default_value("output/"),
             "Destination directory. File will be named based on input file name");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    extract_bubble_chains_from_gfa(
            gfa_path,
            bubble_path,
            assembly_summary_path,
            output_dir);

    return 0;
}
