/* 
 * File:   extractor.h
 * Author: vitor
 *
 * Created on March 31, 2018, 11:45 AM
 */

#include <string>
#include <vector>
#include <iostream>

using namespace std;

enum cartridge_type {EXTRACTION, INDEXING};
enum extension_type {CSV, PROGRAM, FASTBIT, POSTGRES_RAW};

class Extractor{
protected:
    string tag;
    string set_tag;
    cartridge_type cartridge;
    extension_type extension;
    vector<Attribute> attributes;
    
public:
    Extractor(string tag, string set_tag, cartridge_type cartridge, extension_type extension){
        this->tag = tag;
        this->set_tag = set_tag;
        this->cartridge = cartridge;
        this->extension = extension;
    }
    
    void add_attribute(string name, attribute_type type);
    void add_attributes(vector<string> names, vector<attribute_type> types);
    
    string get_tag();
    string get_set_tag();
    Attribute& get_attribute_by_name(string name);
    string get_specification();
};