#include "VCFReader.hpp"
#include <iostream>

using std::cout;

int main() {
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_vcf_path = "/data/test.vcf";
    path absolute_vcf_path = project_directory / relative_vcf_path;

    VCFReader reader(absolute_vcf_path);

    map <string, vector <Variant> > variants;
    reader.read_all(variants);

    for (auto& [chromosome, chromosome_variants]: variants) {
        for (auto& v: chromosome_variants) {
            cout << v.to_string() << '\n';
        }
    }

    return 0;
}

