/* 
 * File:   raw_data_indexer.h
 * Author: vitor
 *
 * Created on April 13, 2018, 8:51 AM
 */

#ifndef RAW_DATA_INDEXER_H
#define RAW_DATA_INDEXER_H

#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>

#include "extractor_enum.h"
#include "attribute_enum.h"

using namespace std;

class RawDataIndexer {
//    command_line << std::getenv("DFANALYZER_DIR")
//                        << "/bin/RDI OPTIMIZED_FASTBIT:INDEX extractor_" << to_string(t_step) 
//                        << " " << string(path) << "/rde/" << to_string(t_step)
//                        << " extractor_"  << to_string(t_step) + ".data"
//                        << " [u:NUMERIC,v:NUMERIC,w:NUMERIC,p:NUMERIC,x:NUMERIC,y:NUMERIC,z:NUMERIC]"
//                        << " -delimiter=\";\" -bin=\"" 
//                        << std::getenv("FASTBIT_DIR") << "/bin\"";        
protected:
    cartridge_type method;
    string cartridge = "INDEX";
    string extractor_tag;
    string path;
    string file_name_with_extracted_data;
    vector<string> attribute_names;
    vector<attribute_type> attribute_types;
    string extra_arguments;

    vector<string> values_of_attribute_types = {"TEXT", "NUMERIC", "RDFILE"};
    vector<string> values_of_cartridge_types = {"CSV", "PROGRAM", "FASTBIT", "OPTIMIZED_FASTBIT", "POSTGRES_RAW"};

public:

    RawDataIndexer(cartridge_type cartridge, string extractor_tag, 
            string path, string file_name_with_extracted_data, 
            vector<string> attribute_names, vector<attribute_type> attribute_types,
            string extra_arguments) {
        this->cartridge = cartridge;
        this->extractor_tag = extractor_tag;
        this->set_attribute_names(attribute_names);
        this->set_attribute_types(attribute_types);
    }

    void set_attribute_names(vector<string> attribute_names);
    void set_attribute_types(vector<attribute_type> attribute_types);

    int run();

    string get_command_line();
    string get_attributes_as_string();
};

#endif /* RAW_DATA_INDEXER_H */

