
#include "GFAReader.hpp"
#include <iostream>
#include <fstream>
#include <string>

using std::stoi;
using std::cout;
using std::runtime_error;


GFAReader::GFAReader(path gfa_path){
    this->gfa_path = gfa_path;

    // Test file
    ifstream test_stream(this->gfa_path);
    if (not test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->gfa_path.string());
    }

}


void GFAReader::index() {
    ifstream gfa_file(this->gfa_path);
    uint64_t line_bytes = 0;
    char c;

    while (gfa_file.get(c)){
        if (c == '\n'){


            line_bytes = 0;
        }
        else{
            if (line_bytes == 0) {
                if (c == 'S'){

                }
                else if (c == 'L'){

                }
            }

            line_bytes++;
        }
    }
}
