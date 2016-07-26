#include <iostream>
#include <string>

#include "provenance_object.h"
#include "../rapidjson/document.h"
#include "../rapidjson/filewritestream.h"
#include "../rapidjson/stringbuffer.h"
#include "../rapidjson/writer.h"
#include "../rapidjson/prettywriter.h"
#include <cstdio>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>

using namespace std;
using namespace rapidjson;

class Task : public ProvenanceObject {
    int subID;
    string dataflow;
    string transformation;
    string workspace;
    map<string, vector<string>> sets; // hash_map<tag,elements<attribute values>>
    string resource;
    enum TaskStatus { RUNNING, FINISHED, FINISHED_WITH_ERROR};
    TaskStatus status;
//    vector<Performance> performances;

public:

    Task(int ID) : ProvenanceObject(ID) {
    };

    void writeJSON(string filename);

    int getSubID() const {
        return subID;
    }

    void setSubID(int subID) {
        this->subID = subID;
    }

};