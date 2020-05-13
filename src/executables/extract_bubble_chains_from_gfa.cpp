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


void write_chain_segments_to_output_gfa(
        BubbleChainComponent& chain_component,
        string_bimap& node_complements,
        GFAReader& gfa_reader,
        string& gfa_line,
        ofstream& output_gfa){

    // Write all the segments to a file
    for (auto& segment: chain_component.segments) {
        try {
            gfa_reader.read_line(gfa_line, gfa_reader.sequence_line_indexes_by_node.at(segment));
        }
        catch (exception& e){
            cerr << "Could not find node in GFA: " << segment << '\n';
            throw e.what();
        }

        output_gfa << gfa_line;
    }
}


void write_all_chains_to_output_gfa(
        vector <vector <BubbleChainComponent> >& chains,
        string_bimap& node_complements,
        GFAReader& gfa_reader,
        ofstream& output_gfa){

    string line;
    for (auto& chain: chains){
        for (auto& component: chain){
            write_chain_segments_to_output_gfa(
                    component,
                    node_complements,
                    gfa_reader,
                    line,
                    output_gfa);
        }
    }
}


void find_single_stranded_chains_from_bubble_chains(
        vector <vector <BubbleChainComponent> >& chains,
        vector <vector <BubbleChainComponent> >& single_stranded_chains,
        unordered_map <string, size_t>& chain_indexes_by_node_ids,
        string_bimap& node_complements,
        unordered_set <string>& single_stranded_nodes){

    size_t i_complement = 0;
    string start_id;
    string start_id_complement;
    unordered_set <size_t> found_chains;
    bool not_found;

    for (size_t i=0; i<chains.size(); i++){
        start_id = chains[i][0].segments[0];
        auto start_id_left = node_complements.left.find(start_id);
        auto start_id_right = node_complements.right.find(start_id);

        // Find the start_id complement id
        if (start_id_left != node_complements.left.end()){
            start_id_complement = start_id_left->get_right();
        }
        else if (start_id_right != node_complements.right.end()){
            start_id_complement = start_id_right->get_left();
        }

        // Use the complement id to search for the complement bubble chain
        i_complement = chain_indexes_by_node_ids.at(start_id_complement);

        // Verify complementary chains are the same size
        if (chains[i].size() != chains[i_complement].size()){
            string error_message;
            error_message += "\tChain IDs (A B): " + std::to_string(chains[i][0].id) + " " + std::to_string(chains[i_complement][0].id) + "\n";
            error_message += "\tNode IDs (A B): " + start_id + " " + start_id_complement + "\n";
            error_message += "\tChain sizes (A B): " + std::to_string(chains[i].size()) + " " + std::to_string(chains[i_complement].size()) + "\n\n";

            for (auto& item: chains[i]){
                error_message += item.to_string();
                error_message += '\n';
            }

            error_message += "\n\n";

            for (auto& item: chains[i_complement]){
                error_message += item.to_string();
                error_message += '\n';
            }

            throw runtime_error("ERROR: bubble chain complement does not match size:\n" + error_message);
        }

        // Check if either of these complementary chains have been found
        not_found = (found_chains.find(i) == found_chains.end()) and (found_chains.find(i_complement) == found_chains.end());

//        cout << "\tChain IDs (A B): " + std::to_string(chains[i][0].id) + " " + std::to_string(chains[i_complement][0].id) + "\n";
//        cout << "\tNode IDs (A B): " + start_id + " " + start_id_complement + "\n";
//        cout << "\tChain sizes (A B): " + std::to_string(chains[i].size()) + " " + std::to_string(chains[i_complement].size()) + "\n";
//        cout << (found_chains.find(i) == found_chains.end()) << " " << (found_chains.find(i_complement) == found_chains.end()) << "\n\n";

        if (not_found){
            single_stranded_chains.push_back(chains[i]);

            // Save the first chain and put it + its complement in a set so they aren't added again
            found_chains.insert(i);
            found_chains.insert(i_complement);

            // Save all the names of the segments so that later the Linkages from the GFA can be found
            for (auto& component: chains[i]){
                for (auto& node: component.segments){
                    single_stranded_nodes.emplace(node);
                }
            }
        }
    }
}


void read_bubble_chains_from_csv(
        ifstream& bubble_chain_file,
        path bubble_path,
        vector <vector <BubbleChainComponent> >& chains,
        unordered_map <string, size_t>& chain_indexes_by_node_ids){

    ///
    /// Parse the bubble chain CSV with the format:
    ///     Chain,Circular,Position,Segment0,Segment1,Segment2,Segment3,Segment4,
    ///
    string line;
    uint64_t l = 0;

    BubbleChainComponent chain_component;
    BubbleChainComponent previous_chain_component;

    string gfa_line;
    string node_name;
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

        // When the chain ends
        if (chain_component.id != previous_chain_component.id and l > 1) {
            if (chain.size() > 1) {
                chains.emplace_back(chain);

                // Record the first and last nodes in this chain, so that its complement can later be identified
                for (auto &s: chain[0].segments) {
                    chain_indexes_by_node_ids.emplace(s, chains.size() - 1);
                }
                for (auto &s: chain[chain.size() - 1].segments) {
                    chain_indexes_by_node_ids.emplace(s, chains.size() - 1);
                }
            }

            chain.resize(0);
        }

        chain.push_back(chain_component);
        previous_chain_component = chain_component;
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

    vector <vector <BubbleChainComponent> > chains;
    vector <vector <BubbleChainComponent> > single_stranded_chains;
    unordered_map <string, size_t> chain_indexes_by_node_ids;
    unordered_set <string> single_stranded_nodes;

    read_bubble_chains_from_csv(
            bubble_chain_file,
            bubble_path,
            chains,
            chain_indexes_by_node_ids);

    find_single_stranded_chains_from_bubble_chains(
            chains,
            single_stranded_chains,
            chain_indexes_by_node_ids,
            node_complements,
            single_stranded_nodes);

    write_all_chains_to_output_gfa(
            single_stranded_chains,
            node_complements,
            gfa_reader,
            output_gfa);

    gfa_reader.write_link_subset_to_file(single_stranded_nodes, output_gfa);
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

