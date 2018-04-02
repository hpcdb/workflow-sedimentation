/* 
 * File:   set.h
 * Author: vitor
 *
 * Created on March 31, 2018, 11:27 AM
 */

#include "attribute.h"
#include "extractor.h"

#include <string>
#include <vector>

using namespace std;

class Set{
protected:
    string tag;
    vector<Attribute> attributes;
    vector<Extractor> extractors;
        
public:
    Set(string tag){
        this->tag = tag;
    }
        
    Attribute& add_attribute(string name, attribute_type type);
    vector<Attribute>& add_attributes(vector<string> names, vector<attribute_type> types);
    Extractor add_extractor(string tag);
    
    Extractor get_extractor_by_tag(string tag);
    string get_tag();
    
    string get_specification();
};