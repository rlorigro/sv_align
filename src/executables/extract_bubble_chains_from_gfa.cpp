#include "BubbleChain.hpp"
#include "GFAReader.hpp"
#include "boost/bimap.hpp"
#include "boost/program_options.hpp"
#include <iostream>
#include <stdexcept>
#include <unordered_set>


using std::cout;
using std::cerr;
using std::getline;
using std::ofstream;
using std::exception;
using std::runtime_error;
using std::unordered_map;
using std::experimental::filesystem::create_directories;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using boost::bimap;

typedef bimap<string,string> string_bimap;
typedef string_bimap::value_type bimap_pair;


void parse_nodes_from_assembly_summary_line(string& line, string_bimap& node_complements){
    ///
    /// Parse a line of the AssemblySummary.csv with the following format:
    ///     Rank,EdgeId,EdgeIdRc,Length,CumulativeLength,LengthFraction,CumulativeFraction
    ///     0,800526,800527,169589,169589,5.80388e-05,5.80388e-05
    ///

    string token;
    uint64_t n_commas = 0;
    string forward_node;
    string reverse_node;

    for (auto& c: line){
        if (c == ','){
            if (n_commas == 1){
                forward_node = token;
            }
            else if (n_commas == 2){
                reverse_node = token;
                node_complements.insert(bimap_pair(forward_node, reverse_node));
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


void extract_node_sets_from_assembly_summary(path assembly_summary_path, string_bimap& node_complements) {
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

        parse_nodes_from_assembly_summary_line(line, node_complements);

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

    string_bimap node_complements;
    extract_node_sets_from_assembly_summary(assembly_summary_path, node_complements);

    GFAReader gfa_reader(gfa_path);
    gfa_reader.map_sequences_by_node();

    string line;
    uint64_t l = 0;
    uint64_t n_forward = 0;
    uint64_t n_reverse = 0;

    BubbleChainComponent chain_component;
    BubbleChainComponent previous_chain_component;

    ///
    /// Parse the bubble chain CSV with the format:
    ///     Chain,Circular,Position,Segment0,Segment1,Segment2,Segment3,Segment4,
    ///
    string gfa_line;
    string node_name;
    unordered_set <string> forward_nodes;
    while (getline(bubble_chain_file, line)){
        // Skip header line
        if (l == 0){
            l++;
            continue;
        }

        // Catch empty lines if they exist?
        if (line.empty()){
            throw runtime_error("ERROR: empty line found in file: " + bubble_path.string() + " at line index " + std::to_string(l));
        }

        // Record the stats about this component in the bubble chain (either a haploid segment or polyploid segment set)
        parse_line_as_bubble_chain_component(line, chain_component);

        // When the chain ends
        if (chain_component.chain != previous_chain_component.chain){
            cout << n_forward << " " << n_reverse << '\n';
            n_forward = 0;
            n_reverse = 0;
        }

        // Write all the segments to a file
        for (auto& segment: chain_component.segments) {
            auto forward_complement = node_complements.left.find(segment);
            auto reverse_complement = node_complements.right.find(segment);

            // If there is a key in the forward complement side of the bimap (the node is already forward)
            if (forward_complement != node_complements.left.end()) {
                n_forward++;
                gfa_reader.read_line(gfa_line, gfa_reader.sequence_line_indexes_by_node.at(segment));
                output_gfa << gfa_line;
            }
            // If there is a key in the reverse complement side of the bimap (the node is reverse, needs to flip)
            else if (reverse_complement != node_complements.right.end()) {
                n_reverse++;
                segment = reverse_complement->get_right();

                try {
                    gfa_reader.read_line(gfa_line, gfa_reader.sequence_line_indexes_by_node.at(segment));
                }
                catch (exception& e){
                    cerr << "Could not find node in GFA: " << segment << '\n';
                    throw e.what();
                }

                output_gfa << gfa_line;
            }
            else{
                throw runtime_error("ERROR: node not in set of known forward/reverse nodes from AssemblySummary.csv: " + segment);
            }
        }

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

