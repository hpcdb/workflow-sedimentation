
#include "set.h"

void Set::add_attribute(string name, attribute_type type) {
    Attribute new_attribute = Attribute(name, type);
    this->attributes.push_back(new_attribute);
}

void Set::add_attributes(vector<string> names, vector<attribute_type> types) {
    if (names.size() == types.size()) {
        for (int index = 0; index < names.size(); index++) {
            string name = names.at(index);
            attribute_type type = types.at(index);
            this->add_attribute(name, type);
        }
    }
}

Extractor& Set::add_extractor(string extractor_tag, cartridge_type cartridge, extension_type extension) {
    Extractor new_extractor = Extractor(extractor_tag, this->tag, cartridge, extension);
    this->extractors.push_back(new_extractor);
}

Extractor& Set::add_extractor(string extractor_tag, cartridge_type cartridge, extension_type extension,
        string attribute_name, attribute_type attribute_type) {
    Extractor& extractor = this->add_extractor(extractor_tag, cartridge, extension);
    extractor.add_attribute(attribute_name, attribute_type);
    this->extractors.push_back(extractor);
}

Extractor& Set::add_extractor(string extractor_tag, cartridge_type cartridge, extension_type extension, 
        vector<string> attribute_names, vector<attribute_type> attribute_types) {
    Extractor& extractor = this->add_extractor(extractor_tag, cartridge, extension);
    extractor.add_attributes(attribute_names, attribute_types);
}

string Set::get_tag() {
    return this->tag;
}

string Set::get_specification() {
    string attribute_names = "";
    string attribute_types = "";

    for (Attribute attribute : this->attributes) {
        if (!attribute_names.empty() and !attribute_types.empty()) {
            attribute_names += ";";
            attribute_types += ";";
        }
        attribute_names += attribute.get_name();
        attribute_types += attribute.get_type();
    }

    attribute_names = "{" + attribute_names + "}";
    attribute_types = "{" + attribute_types + "}";

    return "dataset(" + this->tag + "," +
            attribute_names + "," + attribute_types + ")";
}