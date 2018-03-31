/* 
 * File:   element.h
 * Author: vitor
 *
 * Created on March 31, 2018, 10:32 AM
 */

#include <string>
#include <vector>

using namespace std;

class Element{
protected:
    vector<string> values;
    
public:    
    void add_value(string value);
    
    vector<string> get_values();
    string get_values_as_string();
    
};