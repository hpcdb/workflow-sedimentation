#include <iostream>
#include <string>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <vector>
#include "attribute.h"

using namespace std;

class Indexer {
    string method = "EXTRACTION";
    string extension = "CSV";
    string name = "";
    string path = "";
    string filename = "";
    string delimeter = ",";
    string bin = "";
    string extraArguments = "";
    vector<Attribute> attributes;
    string commandLine = "";

public:

    Indexer(string commandLine, string name) {
        this->name = name;
        this->commandLine = commandLine;
    };

    Indexer(string commandLine, string method, string extension, string name) {
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

    void setBin(string dir){
        this->bin = dir;
    }

    void setExtraArguments(string arguments){
        this->extraArguments = arguments;
    }

    void index(string path, string filename, int indexerID) {
        char buffer[2048];

        if(extension.compare("FASTBIT") == 0 || extension.compare("OPTIMIZED_FASTBIT") == 0){
            sprintf(buffer, "mkdir index/%d; cd index/%d; mv ../../%s .; cp -rf ../../FastBit .;", indexerID, indexerID, filename.c_str());
            sprintf(buffer, "%s%s%s:INDEX %s %s/index/%d %s [", buffer, commandLine.c_str(), extension.c_str(), name.c_str(), path.c_str(), indexerID, filename.c_str());
        }else{
            sprintf(buffer, "%s%s:INDEX %s %s %s [", commandLine.c_str(), extension.c_str(), name.c_str(), path.c_str(), filename.c_str());
        }
        
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
        sprintf(buffer, "%s -performance=%s/%s", buffer, path.c_str(), "prov/indexing/rdi-component.prov");
        if(extension.compare("FASTBIT") == 0 || extension.compare("OPTIMIZED_FASTBIT") == 0){
            sprintf(buffer, "%s -bin=\"%s\"", buffer, this->bin.c_str());
            sprintf(buffer, "%s -option=%s", buffer, this->extraArguments.c_str());
            sprintf(buffer, "%s;cd ../..", buffer);
        }
	    cout << buffer << endl;
        system(buffer);
    }

};