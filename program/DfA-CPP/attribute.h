/* 
 * File:   attribute.h
 * Author: vitor
 *
 * Created on March 31, 2018, 11:30 AM
 */

#include<string>

using namespace std;

enum attribute_type {TEXT, NUMERIC, RDFILE};

class Attribute{
protected:
    string name;
    attribute_type type;
    
public:
    Attribute(string name, attribute_type type){
        this->name = name;
        this->type = type;
    }
    
    string get_name();
    string get_type();
};