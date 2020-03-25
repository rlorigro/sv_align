
#include "VCFReader.hpp"
#include <iostream>
#include <stdexcept>

using std::stoi;
using std::cout;
using std::runtime_error;


string Variant::to_string(){
    return this->chromosome + "\t"
            + std::to_string(this->reference_start) + "\t"
            + std::to_string(this->quality) + "\t"
            + std::to_string(this->pass) + "\t"
            + std::to_string(this->genotype.first) + "/" + std::to_string(this->genotype.second) + "\t"
            + this->alleles[this->genotype.first] + "\t"
            + this->alleles[this->genotype.second];
}


VCFReader::VCFReader(path vcf_path){
    this->vcf_path = vcf_path;

    // Test file
    ifstream test_stream(this->vcf_path);
    if (not test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->vcf_path.string());
    }
}


void VCFReader::parse_genotype(ifstream& vcf_file, char& c, Variant& variant){
    string token;

    while (vcf_file.get(c)){
        if (c == '/' or c == '\n'){
            if (token == "."){
                if (c == '/') {
                    variant.genotype.first = 0;
                    token.resize(0);
                }
                else{
                    variant.genotype.second = 0;
                    break;
                }
            }
            else{
                if (c == '/') {
                    variant.genotype.first = stoi(token);
                    token.resize(0);
                }
                else{
                    variant.genotype.second = stoi(token);
                    break;
                }
            }
        }
        else{
            token += c;
        }
    }
}


void VCFReader::read_all(map <string, vector <Variant> >& variants){
    ///
    /// Parse a vcf with the format:
    /// #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG005733
    ///

    ifstream vcf_file(this->vcf_path);
    string token;
    char c;

    Variant variant;
    uint8_t n_separators = 0;

    while(vcf_file.get(c)){
        // If this line is a header skip until newline character
        if (c == '#'){
            while(c != '\n'){
                vcf_file.get(c);
            }
            token.resize(0);
            continue;
        }

        // If this line contains data, parse it
        if ((c == '\t') or (c == '\n')){
            if (n_separators == 0){
                variant.chromosome = token;
            }
            else if (n_separators == 1){
                variant.reference_start = uint64_t(stoi(token));
            }
            else if (n_separators == 3 or n_separators == 4){
                variant.alleles.push_back(token);
            }
            else if (n_separators == 5){
                variant.quality = uint8_t(stoi(token));
            }
            else if (n_separators == 6){
                if (token == "PASS"){
                    variant.pass = true;
                }
                else{
                    variant.pass = false;
                }
            }
            else if (n_separators == 8){
                parse_genotype(vcf_file, c, variant);
            }

            n_separators++;
            token.resize(0);

            // Reset variant at end of line
            if (c == '\n'){
                // Add a new vector for this chromosome if not in map
                if (variants.count(variant.chromosome) == 0){
                    variants.emplace(variant.chromosome, vector<Variant>());
                }

                variants.at(variant.chromosome).push_back(variant);
                variant = {};
                n_separators = 0;
            }
        }
        else {
            token += c;
        }
    }
}
