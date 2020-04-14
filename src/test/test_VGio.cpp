#include "vg/io/protobuf_iterator.hpp"
#include "vg/vg.pb.h"
#include <experimental/filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <sstream>

using std::experimental::filesystem::path;
using std::ifstream;
using std::string;
using std::pair;
using std::cerr;
using std::stringstream;
using vg::Alignment;


/// Deconstruct a virtual offset into its component parts
static pair<size_t, size_t> unvo(int64_t virtual_offset) {
    pair<size_t, size_t> to_return;

    to_return.first = virtual_offset >> 16;
    to_return.second = virtual_offset & 0xFFFF;
    return to_return;
}


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_gam_path = "/data/HG00733_chr14_85915451_haps_vs_shasta_asm.gam";
    path absolute_gam_path = project_directory / relative_gam_path;

    ifstream datastream(absolute_gam_path);

    for (vg::io::ProtobufIterator<Alignment> it(datastream); it.has_current(); it.advance()) {
        auto vo_parts = unvo(it.tell_group());

        cerr << "Found item " << " at VO " << it.tell_group()
             << " = " << vo_parts.first << ", " << vo_parts.second << '\n';
    }

    return 0;
}
