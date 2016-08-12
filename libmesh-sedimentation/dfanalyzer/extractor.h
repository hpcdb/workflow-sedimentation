#include <iostream>
#include <string>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <vector>
#include "attribute.h"

using namespace std;

class Extractor {
    string method = "EXTRACTION";
    string extension = "CSV";
    string name = "";
    string path = "";
    string filename = "";
    string delimeter = ",";
    vector<Attribute> attributes;
    string commandLine = "";

public:

    Extractor(string commandLine, string name) {
        this->name = name;
        this->commandLine = commandLine;
    };

    Extractor(string commandLine, string method, string extension, string name) {
        this->method = method;
        this->extension = extension;
        this->name = name;
        this->commandLine = commandLine;
    };

    void addAttribute(string name, string type, bool key) {
        Attribute att(name, type, key);
        this->attributes.push_back(att);
    }
    
    void setDelimiter(string delimeter){
        this->delimeter = delimeter;
    }

    void extract(string path, string filename) {
        string execCmd = commandLine + " " + method + ":" + extension + ":EXTRACT " + name + " " + path + " " + filename + " [";
        
        bool first = true;
        for (Attribute att : this->attributes) {
            if(first){
                first = false;
            }else{
                execCmd += ",";
            }
            execCmd += att.GetName() + ":" + att.GetType();
            if(att.IsKey()){
                execCmd += ":key";
            }
        }
        execCmd += "] -delimiter=\"" + this->delimeter + "\"";
        system(strdup(execCmd.c_str()));
    }

};