#include "VCFReader.hpp"
#include "FastaReaderLite.hpp"
#include "boost/program_options.hpp"
#include <iostream>

using std::cout;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void generate_haploblocks_from_vcf(path ref_fasta_path, path vcf_path, path output_dir){
    FastaReaderLite fasta_reader(ref_fasta_path);

    vector <pair <string,string> > sequences;
    fasta_reader.read_all(sequences);

    VCFReader vcf_reader(vcf_path);

    map <string, vector <Variant> > variants;
    vcf_reader.read_all(variants);

    for (auto& [chromosome, chromosome_variants]: variants) {
        for (auto& v: chromosome_variants) {
            
        }
    }
}


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path vcf_path;
    path output_dir;
    uint32_t flank_size;

    options_description options("Arguments");

    options.add_options()
            ("ref",
             value<path>(&ref_fasta_path),
             "File path of reference FASTA file containing REFERENCE sequences")

            ("vcf",
             value<path>(&vcf_path),
             "File path of reference VCF file containing variants to be converted to haploblocks")

            ("output_dir",
             value<path>(&output_dir)->
             default_value("output/"),
             "Destination directory. File will be named based on input file name")

            ("flank_size",
             value<uint32_t>(&flank_size)->
             default_value(100),
             "Destination directory. File will be named based on input file name");

    // Store options in a map and apply values to each corresponding variable
    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);
    notify(vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cout << options << "\n";
        return 0;
    }

    generate_haploblocks_from_vcf(
            ref_fasta_path,
            vcf_path,
            output_dir);

    return 0;
}

