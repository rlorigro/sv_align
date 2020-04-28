#ifndef SV_ALIGN_GFAREADER_H
#define SV_ALIGN_GFAREADER_H

#include <experimental/filesystem>
#include <fstream>
#include <string>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>


using std::experimental::filesystem::path;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::pair;
using std::set;
using std::unordered_set;
using std::map;
using std::unordered_map;

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
    unordered_map <string, size_t> sequence_line_indexes_by_node;
    unordered_map <string, set <size_t> > link_line_indexes_by_node;
    static const char EOF_CODE;

    /// Methods ///
    GFAReader(path gfa_path);
    ~GFAReader();
    void index();
    void read_index();
    void write_index_to_binary_file();
    void map_sequences_by_node();
    void map_links_by_node();
    void read_line(string& s, size_t index);
    void write_link_subset_to_file(unordered_set<string>& node_subset, ofstream& output_file);
};


#endif //SV_ALIGN_GFAREADER_H
