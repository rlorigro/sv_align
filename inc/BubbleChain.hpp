#ifndef SV_ALIGN_BUBBLECHAIN_HPP
#define SV_ALIGN_BUBBLECHAIN_HPP

#include <string>
#include <vector>
#include <stdexcept>

using std::string;
using std::vector;
using std::runtime_error;


class BubbleChainComponent{
public:
    /// Attributes ///
    uint64_t id;
    bool circular;
    uint64_t position;
    vector <string> segments;

    BubbleChainComponent();
    string to_string();
};

void parse_line_as_bubble_chain_component(string& line, BubbleChainComponent& chain_component);


#endif //SV_ALIGN_BUBBLECHAIN_HPP
