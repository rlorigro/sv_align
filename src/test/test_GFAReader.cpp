#include <GFAReader.hpp>
#include <iostream>

using std::cout;

int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_gfa_path = "/data/test_gfa1.gfa";
    path absolute_gfa_path = project_directory / relative_gfa_path;

    GFAReader reader(absolute_gfa_path);

    for (auto& [gfa_type_code, indexes]: reader.line_indexes_by_type){
        for (auto& i: indexes) {
            cout << gfa_type_code << '\t' << i << '\t' << reader.line_offsets[i].offset << '\n';
        }
    }

    cout << reader.line_offsets.back().type << '\t' << reader.line_offsets.back().offset << '\n';

    reader.map_sequences_by_node();

    cout << reader.sequence_line_indexes_by_node.at("11") << '\t' << reader.line_offsets[reader.sequence_line_indexes_by_node.at("11")].offset <<'\n';
    cout << reader.sequence_line_indexes_by_node.at("12") << '\t' << reader.line_offsets[reader.sequence_line_indexes_by_node.at("12")].offset <<'\n';
    cout << reader.sequence_line_indexes_by_node.at("13") << '\t' << reader.line_offsets[reader.sequence_line_indexes_by_node.at("13")].offset <<'\n';

    string s;

    reader.read_line(s, reader.sequence_line_indexes_by_node.at("11"));
    cout << s;
    reader.read_line(s, reader.sequence_line_indexes_by_node.at("12"));
    cout << s;
    reader.read_line(s, reader.sequence_line_indexes_by_node.at("13"));
    cout << s;

    return 0;
}


