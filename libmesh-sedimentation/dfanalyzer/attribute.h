/* 
 * File:   attribute.h
 * Author: vitor
 *
 * Created on August 12, 2016, 3:31 PM
 */
#include <string>

using namespace std;

class Attribute {
    string name;
    string type;
    bool key;
    
public:

    Attribute(string name, string type, bool key) {
        this->name = name;
        this->type = type;
        this->key = key;
    };
    
    bool IsKey() const {
        return key;
    }

    string GetName() const {
        return name;
    }

    string GetType() const {
        return type;
    }


};
