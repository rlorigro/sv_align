
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
        cerr << "No index found, generating .gfai for " << this->gfa_path << "\n";

        this->index();
    }
    // If index is found, load it
    else{
        cerr << "Found index, loading from disk: " << this->gfa_index_path << "\n";

        this->read_index();
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

        this->line_indexes.emplace_back(c,offset);
        this->line_indexes_by_type[c].emplace_back(this->line_indexes.size() - 1);
    }

    ::close(file_descriptor);
}


void GFAReader::write_index_to_binary_file(){
    ofstream index_file(this->gfa_index_path);

    for (auto& item: this->line_indexes){
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
                this->line_indexes.emplace_back(gfa_type_code, byte_offset);
                this->line_indexes_by_type[gfa_type_code].emplace_back(this->line_indexes.size() - 1);
                newline = false;
            }
        }
        byte_offset++;
    }

    // Append a placeholder to tell the total length of the file
    this->line_indexes.emplace_back(this->EOF_CODE, byte_offset);
    this->write_index_to_binary_file();
}
