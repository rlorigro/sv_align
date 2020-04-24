#ifndef SV_ALIGN_GFAREADER_H
#define SV_ALIGN_GFAREADER_H

#include <experimental/filesystem>
#include <fstream>
#include <string>
#include <map>


using std::experimental::filesystem::path;
using std::ifstream;
using std::string;
using std::vector;
using std::pair;
using std::map;

class GFAIndex{
public:
    /// Attributes ///
    char type;
    uint64_t offset;

    /// Methods ///
    GFAIndex(char type, uint64_t offset);
};


class GFAReader {
public:
    /// Attributes ///
    path gfa_path;
    path gfa_index_path;
    int gfa_file_descriptor;
    vector <GFAIndex> line_offsets;
    map <char, vector <size_t> > line_indexes_by_type;
    map <string, size_t> sequence_line_indexes_by_node;
    static const char EOF_CODE;

    /// Methods ///
    GFAReader(path gfa_path);
    ~GFAReader();
    void index();
    void read_index();
    void write_index_to_binary_file();
    void map_sequences_by_node();
    void read_line(string& s, size_t index);
};


#endif //SV_ALIGN_GFAREADER_H
