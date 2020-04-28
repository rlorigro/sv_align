#include <GFAReader.hpp>
#include <iostream>
#include <sstream>

using std::cerr;
using std::cerr;
using std::stringstream;

int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_gfa_path = "/data/test_gfa1.gfa";
    path absolute_gfa_path = project_directory / relative_gfa_path;

    GFAReader reader(absolute_gfa_path);

    for (auto& [gfa_type_code, indexes]: reader.line_indexes_by_type){
        for (auto& i: indexes) {
            cerr << gfa_type_code << '\t' << i << '\t' << reader.line_offsets[i].offset << '\n';
        }
    }

    cerr << reader.line_offsets.back().type << '\t' << reader.line_offsets.back().offset << '\n';

    reader.map_sequences_by_node();

    cerr << reader.sequence_line_indexes_by_node.size() << '\n';
    for (auto& item: reader.sequence_line_indexes_by_node){
        cerr << item.first << " " << item.second << '\n';
    }

    cerr << reader.sequence_line_indexes_by_node.at("11") << '\t' << reader.line_offsets[reader.sequence_line_indexes_by_node.at("11")].offset <<'\n';
    cerr << reader.sequence_line_indexes_by_node.at("12") << '\t' << reader.line_offsets[reader.sequence_line_indexes_by_node.at("12")].offset <<'\n';
    cerr << reader.sequence_line_indexes_by_node.at("13") << '\t' << reader.line_offsets[reader.sequence_line_indexes_by_node.at("13")].offset <<'\n';

    string s;
    
    cerr << "TESTING retrieval\n";
    reader.read_line(s, reader.sequence_line_indexes_by_node.at("11"));
    cerr << s;
    reader.read_line(s, reader.sequence_line_indexes_by_node.at("12"));
    cerr << s;
    reader.read_line(s, reader.sequence_line_indexes_by_node.at("13"));
    cerr << s;

    reader.map_links_by_node();

    cerr << reader.link_line_indexes_by_node.size() << '\n';
    for (auto& item: reader.link_line_indexes_by_node){
        for (auto& s: item.second) {
            cerr << item.first << " " << s << '\n';
        }
    }

    unordered_set<string> nodes = {"11","12"};

    ofstream o("test_GFAReader.test");

    reader.write_link_subset_to_file(nodes, o);



    return 0;
}


