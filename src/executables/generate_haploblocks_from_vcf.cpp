#include "VCFReader.hpp"
#include "FastaReaderLite.hpp"
#include "boost/program_options.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::min;
using std::max;
using std::ofstream;
using std::runtime_error;
using std::to_string;
using std::experimental::filesystem::create_directories;
using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;


void generate_haploblocks_from_vcf(path ref_fasta_path, path vcf_path, uint16_t sample_number, uint32_t flank_size, path output_dir){
    path output_filename = "haploblocks.fasta";
    path output_fasta_path = absolute(output_dir) / output_filename;
    create_directories(output_dir);

    ofstream output_fasta(output_fasta_path);

    if (not output_fasta.is_open()){
        throw runtime_error("ERROR: could not create output file: " + output_fasta_path.string());
    }
    cerr << "Writing to " << output_fasta_path << '\n';

    cerr << "Reading Fasta...\n";
    FastaReaderLite fasta_reader(ref_fasta_path);

    vector <pair <string,string> > sequences;
    fasta_reader.read_all(sequences);

    cerr << "Reading VCF...\n";
    VCFReader vcf_reader(vcf_path);

    map <string, vector <Variant> > variants;
    vcf_reader.read_all(variants, sample_number);

    int64_t left_flank_start = 0;
    int64_t right_flank_start = 0;
    int64_t right_flank_size = 0;

    cerr << "Generating Haploblocks...\n";
    for (auto& [chromosome_name, sequence]: sequences) {
        if (variants.count(chromosome_name) == 0){
            cout << "Skipping " << chromosome_name << '\n';
            continue;
        }

        for (auto& variant: variants.at(chromosome_name)) {
            // Haplotype 0
            left_flank_start = variant.reference_start - flank_size - 1;
            left_flank_start = max(int64_t(0), left_flank_start);
            right_flank_start = variant.reference_start - 1 + variant.alleles[0].size();
            right_flank_start = min(int64_t(sequence.size()), right_flank_start);
            right_flank_size = min(int64_t(sequence.size() - right_flank_start), int64_t(flank_size));

            output_fasta << '>' << variant.chromosome << '_' << to_string(variant.reference_start) << "_h0_" << variant.alleles[variant.genotype.first].size() << '\n';
            output_fasta << sequence.substr(left_flank_start,variant.reference_start - left_flank_start - 1)
                         << variant.alleles[variant.genotype.first]
                         << sequence.substr(right_flank_start,right_flank_size) << '\n';

            // Haplotype 1
            left_flank_start = variant.reference_start - flank_size - 1;
            left_flank_start = max(int64_t(0), left_flank_start);
            right_flank_start = variant.reference_start - 1 + variant.alleles[0].size();
            right_flank_start = min(int64_t(sequence.size() - 1), right_flank_start);
            right_flank_size = min(int64_t(sequence.size() - right_flank_start), int64_t(flank_size));

            output_fasta << '>' << variant.chromosome << '_' << to_string(variant.reference_start) << "_h1_" << variant.alleles[variant.genotype.second].size() << '\n';
            output_fasta << sequence.substr(left_flank_start,variant.reference_start - left_flank_start - 1);
            output_fasta << variant.alleles[variant.genotype.second];
            output_fasta << sequence.substr(right_flank_start,right_flank_size) << '\n';

        }
    }
}


int main(int argc, char* argv[]){
    path ref_fasta_path;
    path vcf_path;
    path output_dir;
    uint32_t flank_size;
    uint16_t sample_number;

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
             "Destination directory. File will be named based on input file name")

            ("sample",
             value<uint16_t>(&sample_number)->
             default_value(0),
             "The number of the sample (in order of appearance) to use for generating haploblocks, STARTING FROM 0");

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
            sample_number,
            flank_size,
            output_dir);

    return 0;
}

