#ifndef SV_ALIGN_VCFREADER_HPP
#define SV_ALIGN_VCFREADER_HPP

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


class Variant{
public:
    /// Attributes ///
    string chromosome;
    uint64_t reference_start;
    uint8_t quality;
    bool pass;
    pair <uint8_t, uint8_t> genotype;
    vector <string> alleles;

    /// Methods ///
    string to_string(char separator='\t');
};

class VCFReader {
public:
    /// Attributes ///
    path vcf_path;

    /// Methods ///
    VCFReader(path vcf_path);
    void read_all(map <string, vector <Variant> >& variants, uint16_t sample_number=0);
    void parse_genotype(ifstream& vcf_file, char& c, Variant& variant);
};


#endif //SV_ALIGN_VCFREADER_HPP
