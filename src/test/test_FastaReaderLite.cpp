#include "FastaReaderLite.hpp"
#include <iostream>

using std::cout;

int main() {
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_fasta_path = "/data/test.fasta";
    path absolute_fasta_path = project_directory / relative_fasta_path;

    FastaReaderLite reader(absolute_fasta_path);

    vector <pair <string,string> > sequences;
    reader.read_all(sequences);

    for (auto& [header, sequence]: sequences) {
        cout << header << '\n' << sequence << '\n';
    }

    return 0;
}

