
#include "set.h"

Attribute& Set::add_attribute(string name, attribute_type type) {
    Attribute new_attribute = Attribute(name, type);
    this->attributes.push_back(new_attribute);
    return this->attributes.at(this->attributes.size() - 1);
}

vector<Attribute>& Set::add_attributes(vector<string> names, vector<attribute_type> types) {
    if (names.size() == types.size()) {
        for (int index = 0; index < names.size(); index++) {
            string name = names.at(index);
            attribute_type type = types.at(index);
            this->add_attribute(name, type);
        }
    }
    return this->attributes;
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