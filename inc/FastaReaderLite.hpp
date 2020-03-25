#ifndef SV_ALIGN_FASTAREADERLITE_H
#define SV_ALIGN_FASTAREADERLITE_H

#include <experimental/filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <map>

using std::experimental::filesystem::path;
using std::ifstream;
using std::string;
using std::vector;
using std::pair;
using std::map;


class FastaReaderLite {
public:
    /// Attributes ///
    path fasta_path;

    /// Methods ///
    FastaReaderLite(path vcf_path);
    void read_all(vector <pair <string,string> >& sequences);
};


#endif //SV_ALIGN_FASTAREADERLITE_H
