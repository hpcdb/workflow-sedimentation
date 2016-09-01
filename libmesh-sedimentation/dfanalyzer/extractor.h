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

    void setDelimiter(string delimeter) {
        this->delimeter = delimeter;
    }

    void extract(string path, string filename) {
        char* buffer = (char*) malloc(512);
        sprintf(buffer, "#!/bin/bash;%s%s:%s:EXTRACT %s %s %s [", commandLine.c_str(), method.c_str(), extension.c_str(), name.c_str(), path.c_str(), filename.c_str());
        
        bool first = true;
        for (Attribute att : this->attributes) {
            if (first) {
                first = false;
            } else {
                sprintf(buffer, "%s,", buffer);
            }
            sprintf(buffer, "%s%s:%s", buffer, att.GetName().c_str(), att.GetType().c_str());
            if (att.IsKey()) {
                sprintf(buffer, "%s:key", buffer);
            }
        }
        sprintf(buffer, "%s] -delimiter=\"%s\"", buffer, this->delimeter.c_str());
        cout << buffer << endl;
        system(buffer);
        free(buffer);
    }

};