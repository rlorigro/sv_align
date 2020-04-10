#ifndef SV_ALIGN_GFAREADER_H
#define SV_ALIGN_GFAREADER_H

#include <experimental/filesystem>
#include <fstream>
#include <string>

using std::experimental::filesystem::path;
using std::ifstream;
using std::string;

class GFAIndex{
public:
    /// Attributes ///
    path gfa_path;

};

class GFAReader {
public:
    /// Attributes ///
    path gfa_path;

    /// Methods ///
    GFAReader(path gfa_path);
    void index();
};


#endif //SV_ALIGN_GFAREADER_H
