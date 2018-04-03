#include "element.h"
#include <boost/algorithm/string/join.hpp>

void Element::add_value(string value){
    this->values.push_back(value);
}

void Element::set_values(vector<string> values){
    this->values.assign(values.begin(), values.end());
}

vector<string>& Element::get_values(){
    return this->values;
}

string Element::get_values_as_string(){
    return boost::algorithm::join(this->values, ";");
    
    stringstream values;
    for(string value : this->values){
        if(!values.str().empty()){
            values << ";";
        }
        values << value;
    }
    return values.str();
}