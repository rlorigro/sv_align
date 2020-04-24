#include "BubbleChain.hpp"


string BubbleChainComponent::to_string(){
    string s;
    s = std::to_string(this->chain) + '\t' + std::to_string(this->circular) + '\t' + std::to_string(this->position);
    for (auto& segment: this->segments){
        s += '\t' + segment;
    }

    return s;
}


BubbleChainComponent::BubbleChainComponent(){
    this->chain = 0;
    this->circular = false;
    this->position = 0;
}


bool parse_string_as_bool(string& s){
    if (s == "Yes"){
        return true;
    }
    else if (s == "No"){
        return false;
    }
    else{
        throw runtime_error("ERROR: cannot parse string as boolean: " + s);
    }
}


void parse_line_as_bubble_chain_component(string& line, BubbleChainComponent& chain_component){
    chain_component = {};
    uint64_t n_commas = 0;
    string token;

    for (auto& c: line){
        if (c == ',') {
            if (n_commas == 0){
                chain_component.chain = stoi(token);
            }
            else if (n_commas == 1){
                chain_component.circular = parse_string_as_bool(token);
            }
            else if (n_commas == 2){
                chain_component.position = stoi(token);
            }
            else if (n_commas >= 3 and not token.empty()){
                chain_component.segments.push_back(token);
            }
            n_commas++;
            token.resize(0);
        }
        else {
            token += c;
        }
    }
}


