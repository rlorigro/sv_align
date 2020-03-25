#include "FastaReaderLite.hpp"
#include <iostream>
#include <stdexcept>

using std::cout;
using std::getline;
using std::runtime_error;


FastaReaderLite::FastaReaderLite(path fasta_path){
    this->fasta_path = fasta_path;

    // Test file
    ifstream test_stream(this->fasta_path);
    if (not test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->fasta_path.string());
    }
}


void FastaReaderLite::read_all(vector <pair <string,string> >& sequences){
    ifstream fasta_file(this->fasta_path);
    string line;
    string header;
    string sequence;

    while(getline(fasta_file,line)){
        if (line[0] == '>'){
            // If not the first iteration, append a sequence element (pair of header + nt) to the vector
            if (not header.empty()){
                sequences.emplace_back(header, sequence);
            }

            header = line.substr(1,line.size()-1);
        }
        else{
            sequence = line.substr(0,line.size()-1);
        }
    }

    // Tag on the last sequence before the end of file
    if (not sequence.empty()){
        sequences.emplace_back(header, sequence);
    }

}