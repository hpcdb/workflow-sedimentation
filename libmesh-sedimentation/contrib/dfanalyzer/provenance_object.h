#include <iostream>
#include <string>

#include "../rapidjson/document.h"
#include "../rapidjson/filewritestream.h"
#include "../rapidjson/stringbuffer.h"
#include "../rapidjson/writer.h"
#include "../rapidjson/prettywriter.h"
#include <cstdio>
#include <sstream>
#include <fstream>

using namespace std;
using namespace rapidjson;

class ProvenanceObject {
protected:
    int ID = 0;

public:

    ProvenanceObject(int newID);

    virtual void sendRequest(string dfa_hostname);
    
    int getID() const {
        return ID;
    }

    void setID(int ID) {
        this->ID = ID;
    }
};