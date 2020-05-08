#include "BubbleChain.hpp"
#include "GFAReader.hpp"
#include "vg/vg.pb.h"
#include "vg/io/protobuf_iterator.hpp"
#include "boost/program_options.hpp"
#include <string>
#include <experimental/filesystem>
#include "boost/bimap.hpp"

using std::ifstream;
using std::string;
using std::pair;
using std::cout;
using std::cerr;
using std::stoi;
using std::to_string;
using std::unordered_set;
using std::runtime_error;
using std::experimental::filesystem::path;
using std::experimental::filesystem::create_directories;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using vg::Alignment;
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
        throw runtime_error("ERROR: could not open assembly summary file: " + assembly_summary_path.string());
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


void extract_haplotype_info_from_read_name(string& read_name, uint64_t& haplotype, uint64_t& length){
    uint64_t n_separators = 0;

    for (int64_t i=read_name.size()-1; i>=0; i--){
        if (read_name[i] == '_'){
            if (n_separators == 0){
                length = stoi(read_name.substr(i+1,read_name.size()-i));
                haplotype = stoi(string(1,read_name[i-1]));

                cerr << read_name << " " << read_name.substr(i+1,read_name.size()-i) << " " << length << " " << haplotype << '\n';
            }
            n_separators++;
        }
    }
}


void measure_sv_sensitivity(path gfa_path, path gam_path, path bubble_path, path assembly_summary_path, path output_dir){
    GFAReader gfa_reader(gfa_path);
    gfa_reader.map_sequences_by_node();
    gfa_reader.map_links_by_node();

    create_directories(output_dir);
    ifstream bubble_chain_file(bubble_path);
    ifstream gam_file(gam_path);

    if (not bubble_chain_file.is_open()){
        throw runtime_error("ERROR: could not open bubble chain file: " + bubble_path.string());
    }

    if (not gam_file.is_open()){
        throw runtime_error("ERROR: could not open bubble chain file: " + bubble_path.string());
    }

    path output_path = output_dir / ("bubble_stats_" + gam_path.filename().string());
    output_path.replace_extension("csv");
    ofstream output_file(output_path);

    string_bimap node_complements;
    extract_node_sets_from_assembly_summary(assembly_summary_path, node_complements);

    string line;
    uint64_t l = 0;

    BubbleChainComponent chain_component;
    BubbleChainComponent previous_chain_component;

    ///
    /// Parse the bubble chain CSV with the format:
    ///     Chain,Circular,Position,Segment0,Segment1,Segment2,Segment3,Segment4,
    ///
    string gfa_line;
    string node_name;
    unordered_set <string> is_bubble;
    vector <BubbleChainComponent> chain;
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

        // Record the data about this component in the bubble chain (either a haploid segment or polyploid segment set)
        parse_line_as_bubble_chain_component(line, chain_component);

        if (chain_component.segments.size() > 1){
            for (auto& segment: chain_component.segments) {
                is_bubble.insert(segment);

                auto forward_complement = node_complements.left.find(segment);
                auto reverse_complement = node_complements.right.find(segment);
                if (forward_complement != node_complements.left.end()) {
                    is_bubble.insert(forward_complement->get_right());
                }
                else if (reverse_complement != node_complements.right.end()) {
                    is_bubble.insert(reverse_complement->get_left());
                }
            }
        }

        for (auto& segment: chain_component.segments){
            if (segment == "464300"){
                cerr << chain_component.to_string() << '\n';
                cerr << chain_component.segments.size() << '\n';
            }
        }

    }

    node_complements.clear();

    string read_name;
    uint64_t haplotype = 0;
    uint64_t haplotype_length = 0;
    uint16_t n_bubbles = 0;
    uint64_t total_bubble_length = 0;
    bool alignment_in_bubble;
    ifstream datastream(gam_path);
    unordered_set<string> nodes_in_alignment;

    for (vg::io::ProtobufIterator<Alignment> it(datastream); it.has_current(); it.advance()) {
        Alignment& alignment = *it;
        read_name = alignment.name();

        nodes_in_alignment.clear();
        node_name.resize(0);

        total_bubble_length = 0;
        n_bubbles = 0;

        cout << "\nName: " << read_name << '\n';

        for (auto& mapping: alignment.path().mapping()){
            node_name = to_string(mapping.position().node_id());
            nodes_in_alignment.insert(node_name);

            alignment_in_bubble = (is_bubble.find(node_name) != is_bubble.end());

            cout << "Segment: " << node_name << '\n';
            cout << "Is bubble: " << alignment_in_bubble << '\n';
            cout << "Is bubble??: " << (is_bubble.count(node_name) > 0) << '\n';
            cout << "Is bubble??: " << (is_bubble.find(node_name) != is_bubble.end()) << '\n';

            if (alignment_in_bubble){
                n_bubbles++;
                total_bubble_length += gfa_reader.get_sequence_length(node_name);
            }
        }

        extract_haplotype_info_from_read_name(read_name, haplotype, haplotype_length);

        if (haplotype_length > 1){
            output_file << read_name << "," << n_bubbles << "," << haplotype_length << "," << total_bubble_length << '\n';
        }

        if (read_name == "chr14_85915451_h0_1") {
            path subgraph_path = output_dir / (read_name + "_subgraph.gfa");
            ofstream subgraph_gfa_file(subgraph_path);
            gfa_reader.write_subgraph_to_file(nodes_in_alignment, subgraph_gfa_file);
        }
    }
}


int main(int argc, char* argv[]){
    path gfa_path;
    path gam_path;
    path bubble_path;
    path assembly_summary_path;
    path output_dir;

    options_description options("Arguments");

    options.add_options()
            ("gfa",
             value<path>(&gfa_path),
             "File path of GFA file containing reference graph which reads were aligned to")

            ("gam",
             value<path>(&gam_path),
             "File path of GAM file containing SVs aligned to assembly GFA")

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

    measure_sv_sensitivity(
            gfa_path,
            gam_path,
            bubble_path,
            assembly_summary_path,
            output_dir);

    return 0;
}

