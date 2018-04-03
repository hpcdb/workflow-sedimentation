/* 
 * File:   dataset.h
 * Author: vitor
 *
 * Created on March 31, 2018, 10:29 AM
 */

#include "element.h"

#include <string>
#include <vector>

using namespace std;

class Dataset{
protected:
    string tag;
    vector<Element> elements;
    
public:
    Dataset(string tag){
        this->tag = tag;
    }
    
    void add_element(Element element);
    
    string get_tag();
    vector<Element>& get_elements();    
};