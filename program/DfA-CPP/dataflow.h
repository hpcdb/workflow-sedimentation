/* 
 * File:   dataflow.h
 * Author: vitor
 *
 * Created on March 31, 2018, 11:03 AM
 */

#include "transformation.h"
#include "set.h"

#include <string>
#include <vector>

using namespace std;

class Dataflow{
protected:
    string tag;
    vector<Transformation> transformations;
    vector<Set> sets;
    
public:
    Dataflow(string tag){
        this->tag = tag;
    }
    
    int begin(); // return a code with the status of this service request
    int end();
    
    void add_transformation(string tag);
    void add_set(string tag);
};