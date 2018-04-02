
#include "set.h"

Attribute Set::add_attribute(string name, attribute_type type) {
    Attribute new_attribute = Attribute(name, type);
    this->attributes.push_back(new_attribute);
    return new_attribute;
}

string Set::get_tag(){
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