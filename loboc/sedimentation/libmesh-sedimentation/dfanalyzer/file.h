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

class File {
    string path;
    string filename;
    
public:

    File(string path, string filename) {
        this->path = path;
        this->filename = filename;
    };

    string GetPath(){
        return path;
    }
    
    string GetFileName(){
        return filename;
    }

};