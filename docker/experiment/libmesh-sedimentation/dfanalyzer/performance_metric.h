#include "../rapidjson/document.h"
#include "../rapidjson/filewritestream.h"
#include "../rapidjson/stringbuffer.h"
#include "../rapidjson/writer.h"
#include "../rapidjson/prettywriter.h"
#include "date.h"

#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>

#include <chrono>
#include <string>
#include <sstream>

using namespace std;
using namespace rapidjson;

class PerformanceMetric {
    string description;
    string method;
    string startTime;
    string endTime;
    
public:

    PerformanceMetric() {
    };

    string GetDescription() const {
        return description;
    }

    void SetDescription(string description) {
        this->description = description;
    }

    string GetEndTime() const {
        return endTime;
    }

    void SetEndTime(string endTime) {
        this->endTime = endTime;
    }

    string GetMethod() const {
        return method;
    }

    void SetMethod(string method) {
        this->method = method;
    }

    string GetStartTime() const {
        return startTime;
    }

    void SetStartTime(string startTime) {
        this->startTime = startTime;
    }
    
    void IdentifyStartTime();
    
    void IdentifyEndTime();

};