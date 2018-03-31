/* 
 * File:   extractor.h
 * Author: vitor
 *
 * Created on March 31, 2018, 11:45 AM
 */

#include <string>
#include <vector>

using namespace std;

class Extractor{
protected:
    string tag;
    enum cartridge_type {EXTRACTION, INDEXING};
    cartridge_type cartridge;
    enum extension_type {CSV, PROGRAM, FASTBIT, POSTGRES_RAW};
    extension_type extension;
    vector<Attribute> attributes;
    
public:
    Extractor(string tag, cartridge_type cartridge, extension_type extension){
        this->tag = tag;
        this->cartridge = cartridge;
        this->extension = extension;
    }
    
    Attribute add_attribute(string name, attribute_type type);
};