#include "dataset.h"

void Dataset::add_element(Element element){
    this->elements.push_back(element);
}
    
string Dataset::get_tag(){
    return this->tag;
}

vector<Element>& Dataset::get_elements(){
    return this->get_elements();
}