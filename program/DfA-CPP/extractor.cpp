
#include "attribute.h"
#include "extractor.h"

#include <exception>
        
void Extractor::add_attribute(string name, attribute_type type) {
    Attribute new_attribute = Attribute(name, type);
    this->attributes.push_back(new_attribute);
}

void Extractor::add_attributes(vector<string> names, vector<attribute_type> types) {
    if (names.size() == types.size()) {
        for (int index = 0; index < names.size(); index++) {
            string name = names.at(index);
            attribute_type type = types.at(index);
            this->add_attribute(name, type);
        }
    }
}

string Extractor::get_tag() {
    return this->tag;
}

string Extractor::get_set_tag() {
    return this->set_tag;
}

Attribute& Extractor::get_attribute_by_name(string name) {
    try {
        for (int index = 0; index< this->attributes.size(); index++) {
            Attribute& attribute = this->attributes.at(index);
            if (attribute.get_name() == name) {
                return attribute;
            }
        }
    } catch (exception e) {
        cout << "Attribute " + name + " was not captured by extractor " + this->get_tag() << endl;
    }
}

string Extractor::get_specification() {
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

    return "extractor(" + this->tag + "," + this->set_tag + "," +
            attribute_names + "," + attribute_types + ")";
}


