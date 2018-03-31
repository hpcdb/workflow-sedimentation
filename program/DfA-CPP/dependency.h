/* 
 * File:   dependency.h
 * Author: vitor
 *
 * Created on March 31, 2018, 11:15 AM
 */

#include <string>
#include <vector>

using namespace std;

class Dependency{
protected:
    vector<string> transformation_tags;
    vector<vector<int>> transformation_ids; 
    
public:
    Dependency(...){
        
    }
    
    void add_transformation_ids(...);
};