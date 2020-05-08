
#include "GFAReader.hpp"
#include "BinaryIO.hpp"
#include <iostream>
#include <fstream>
#include <string>

using std::stoi;
using std::cout;
using std::ofstream;
using std::runtime_error;


const char GFAReader::EOF_CODE = 'X';

GFAIndex::GFAIndex(char type, uint64_t offset){
    this->type = type;
    this->offset = offset;
}


GFAReader::GFAReader(path gfa_path){
    this->gfa_path = gfa_path;
    this->gfa_index_path = gfa_path;
    this->gfa_index_path.replace_extension("gfai");
    this->gfa_file_descriptor = -1;

    // Test file
    ifstream test_stream(this->gfa_path);
    if (not test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->gfa_path.string());
    }

    // Check if index exists, and generate one if necessary
    if (!exists(this->gfa_index_path)) {
        cerr << "No index found, generating .gfai for " << this->gfa_path << " ... ";

        this->index();
        cerr << "done\n";
    }
    // If index is found, load it
    else{
        cerr << "Found index, loading from disk: " << this->gfa_index_path << " ... ";

        this->read_index();
        cerr << "done\n";
    }
}


GFAReader::~GFAReader(){
    if (fcntl(this->gfa_file_descriptor, F_GETFD)){
        ::close(this->gfa_file_descriptor);
    }
}


void GFAReader::read_index(){
    // Open the input file.
    int file_descriptor = ::open(this->gfa_index_path.c_str(), O_RDONLY);

    // Verify it is working
    if(file_descriptor == -1) {
        throw runtime_error("ERROR: could not read " + this->gfa_index_path.string());
    }

    // Find file size in bytes, calculate number of entries in file
    off_t file_length = lseek(file_descriptor, 0, SEEK_END);
    off_t n_entries = file_length / (sizeof(uint64_t) + sizeof(char));

    // Initialize index to track position in file (AKA cursor)
    off_t byte_index = 0;

    // Read all the chunks and add them to the index data structures
    for (off_t i = 0; i < n_entries; i++){
        char c;
        uint64_t offset;

        pread_value_from_binary(file_descriptor,  c, byte_index);
        pread_value_from_binary(file_descriptor,  offset, byte_index);

        this->line_offsets.emplace_back(c, offset);
        this->line_indexes_by_type[c].emplace_back(this->line_offsets.size() - 1);
    }

    ::close(file_descriptor);
}


void GFAReader::write_index_to_binary_file(){
    ofstream index_file(this->gfa_index_path);

    for (auto& item: this->line_offsets){
        // Write a pair, denoting the line type and the byte offset for the start of that line
        write_value_to_binary(index_file, item.type);
        write_value_to_binary(index_file, item.offset);
    }
}


void GFAReader::index() {
    ifstream gfa_file(this->gfa_path);
    uint64_t byte_offset = 0;
    bool newline = true;
    char gfa_type_code;
    char c;

    // Find all the newlines in the GFA and store each line's byte offset in a vector. Additionally build a map which
    // lists all the positions in the index vector for each line type (e.g. S,L,H,U, etc.), so they can be iterated even
    // if they are not grouped or in order (which is not required by the GFA format spec)
    while (gfa_file.get(c)){
        if (c == '\n'){
            newline = true;
        }
        else{
            if (newline) {
                gfa_type_code = c;
                this->line_offsets.emplace_back(gfa_type_code, byte_offset);
                this->line_indexes_by_type[gfa_type_code].emplace_back(this->line_offsets.size() - 1);
                newline = false;
            }
        }
        byte_offset++;
    }

    // Append a placeholder to tell the total length of the file
    this->line_offsets.emplace_back(this->EOF_CODE, byte_offset);
    this->write_index_to_binary_file();
}


void GFAReader::read_line(string& s, size_t index){
    if (this->gfa_file_descriptor == -1){
        this->gfa_file_descriptor = ::open(this->gfa_path.c_str(), O_RDONLY);
    }

    off_t offset_start = this->line_offsets[index].offset;
    off_t offset_stop = this->line_offsets[index+1].offset;
    off_t length = offset_stop - offset_start;

    pread_string_from_binary(this->gfa_file_descriptor, s, length, offset_start);
}


void GFAReader::map_sequences_by_node(){
    cerr << "Mapping GFA S lines to node names... ";

    ifstream gfa_file(this->gfa_path);
    string token;
    char c = 0;

    // For every sequence line that has been indexed, jump to their offset in the file and read just the node names
    for (size_t i=0; i<this->line_indexes_by_type.at('S').size(); i++){
        auto line_index = this->line_indexes_by_type.at('S')[i];
        auto offset_start = this->line_offsets[line_index].offset + 2;  // Skip the line type character and tab
        token.resize(0);

        gfa_file.seekg(offset_start);

        // Parse the node name
        while (gfa_file.get(c)){
            if (c == '\t'){
                break;
            }
            token += c;
        }

        // Create the mapping that tells where in the file to find each node's sequence data
        this->sequence_line_indexes_by_node[token] = line_index;
        c = 0;
    }

    cerr << "done\n";
}


void GFAReader::map_links_by_node(){
    cerr << "Mapping GFA L lines to node names... ";

    ifstream gfa_file(this->gfa_path);
    string token;
    uint64_t n_separators = 0;
    char c = 0;

    // For every link line that has been indexed, jump to their offset in the file and read just the node names
    for (size_t i=0; i<this->line_indexes_by_type.at('L').size(); i++){
        auto line_index = this->line_indexes_by_type.at('L')[i];
        auto offset_start = this->line_offsets[line_index].offset + 2;  // Skip the line type character and tab
        token.resize(0);
        n_separators = 0;

        gfa_file.seekg(offset_start);

        // Parse the node name
        while (gfa_file.get(c)){
            if (c == '\t'){
                if (n_separators == 0){
                    // Create the mapping that tells where in the file to find each node's link data
                    this->link_line_indexes_by_node[token].insert(line_index);
                }
                else if (n_separators == 2){
                    // Create the mapping that tells where in the file to find each node's link data
                    this->link_line_indexes_by_node[token].insert(line_index);
                    break;
                }
                token.resize(0);
                n_separators++;
            }
            else {
                token += c;
            }
        }

        c = 0;
    }

    cerr << "done\n";
}


void GFAReader::write_link_subset_to_file(unordered_set<string>& node_subset, ofstream& output_file){
    cerr << "Writing GFA L lines to file... ";

    ifstream gfa_file(this->gfa_path);
    uint64_t n_separators = 0;
    bool found_a = false;
    bool found_b = false;
    string token;
    string line;
    char c = 0;

    // For every link line that has been indexed, jump to their offset in the file and read just the node names
    for (size_t i=0; i<this->line_indexes_by_type.at('L').size(); i++){
        auto line_index = this->line_indexes_by_type.at('L')[i];
        auto offset_start = this->line_offsets[line_index].offset;  // Skip the line type character and tab
        token.resize(0);
        n_separators = 0;
        found_a = false;
        found_b = false;

        gfa_file.seekg(offset_start);

        // Parse the node name
        while (gfa_file.get(c)){
            if (c == '\t'){
                if (n_separators == 1){
                    found_a = (node_subset.count(token) != 0);
                }
                else if (n_separators == 3){
                    found_b = (node_subset.count(token) != 0);
                }
                token.resize(0);
                n_separators++;
            }
            else {
                token += c;
            }
            line += c;

            if (c == '\n'){
                break;
            }
        }

        if (found_a and found_b){
            output_file << line;
        }

        line.resize(0);
        c = 0;
    }

    cerr << "done\n";
}


void GFAReader::write_subgraph_to_file(unordered_set <string>& nodes, ofstream& output_gfa){
    string gfa_line;
    for (auto& node_name: nodes){
        this->read_line(gfa_line, this->sequence_line_indexes_by_node.at(node_name));
        output_gfa << gfa_line;
    }

    this->write_link_subset_to_file(nodes, output_gfa);
}


uint64_t GFAReader::get_sequence_length(string node_name){
    auto vector_index = this->sequence_line_indexes_by_node.at(node_name);
    auto start_index = this->line_offsets[vector_index].offset;

    char c;
    uint64_t length = 0;
    uint64_t n_separators = 0;
    ifstream file(this->gfa_path);

    file.seekg(start_index);

    // Count only characters in the sequence (not before or after)
    while(file.get(c)) {
        if (c == '\t' or c == '\n') {
            n_separators++;
        } else {
            if (n_separators == 2) {
                length++;
            }
            else if (n_separators > 2){
                break;
            }
        }
    }

    return length;
}

