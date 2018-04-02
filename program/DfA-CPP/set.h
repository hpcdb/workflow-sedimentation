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

class Set {
protected:
    string tag;
    vector<Attribute> attributes;
    vector<Extractor> extractors;

public:

    Set(string tag) {
        this->tag = tag;
    }

    void add_attribute(string name, attribute_type type);
    void add_attributes(vector<string> names, vector<attribute_type> types);

    Extractor& add_extractor(string extractor_tag, cartridge_type cartridge, extension_type extension);
    Extractor& add_extractor(string extractor_tag, cartridge_type cartridge, extension_type extension,
            string attribute_name, attribute_type attribute_type);
    Extractor& add_extractor(string extractor_tag, cartridge_type cartridge, extension_type extension,
            vector<string> attribute_names, vector<attribute_type> attribute_types);

    string get_tag();
    string get_specification();
};