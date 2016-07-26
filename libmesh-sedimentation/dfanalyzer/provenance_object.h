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
    int ID;

public:

    ProvenanceObject(int newID);

    virtual void writeJSON(string filename);

};