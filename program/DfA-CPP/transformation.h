/* 
 * File:   transformation.h
 * Author: vitor
 *
 * Created on March 31, 2018, 11:07 AM
 */

#include <string>

using namespace std;

class Transformation{
protected:
    string tag;
    
public:
    Transformation(string tag){
        this->tag = tag;
    }
};