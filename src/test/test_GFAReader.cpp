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
            cout << gfa_type_code << '\t' << i << '\t' << reader.line_indexes[i].offset << '\n';
        }
    }

    cout << reader.line_indexes.back().type << '\t' << reader.line_indexes.back().offset << '\n';

    return 0;
}


