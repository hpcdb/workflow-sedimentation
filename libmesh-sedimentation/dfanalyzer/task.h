#include <iostream>
#include <string>

#include "provenance_object.h"
#include "performance_metric.h"
#include "../rapidjson/document.h"
#include "../rapidjson/filewritestream.h"
#include "../rapidjson/stringbuffer.h"
#include "../rapidjson/writer.h"
#include "../rapidjson/prettywriter.h"
#include "../rapidjson/rapidjson.h"
#include <cstdio>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>

using namespace std;
using namespace rapidjson;

class Task : public ProvenanceObject {
protected:
    int subID;
    string dataflow;
    string transformation;
    string workspace;
    map<string, vector<string>> sets; // hash_map<tag,elements<attribute values>>
    string resource;
    string status;
    vector<PerformanceMetric> performances;

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
    
    string getDataflow() const {
        return dataflow;
    }

    void setDataflow(string dataflow) {
        this->dataflow = dataflow;
    }

    vector<PerformanceMetric> getPerformances() const {
        return performances;
    }
    
    void addPerformance(PerformanceMetric perf) {
        this->performances.push_back(perf);
    }

    void setPerformances(vector<PerformanceMetric> performances) {
        this->performances = performances;
    }

    string getResource() const {
        return resource;
    }

    void setResource(string resource) {
        this->resource = resource;
    }

    map<string, vector<string> > getSets() const {
        return sets;
    }
    
    void addSet(string set, vector<string> elements) {
        this->sets.insert(pair<string,vector<string>>(set,elements));
    }

    void setSets(map<string, vector<string> > sets) {
        this->sets = sets;
    }

    string getStatus() const {
        return status;
    }

    void setStatus(string status) {
        this->status = status;
    }

    string getTransformation() const {
        return transformation;
    }

    void setTransformation(string transformation) {
        this->transformation = transformation;
    }

    string getWorkspace() const {
        return workspace;
    }

    void setWorkspace(string workspace) {
        this->workspace = workspace;
    }


};