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
    vector<Element> get_elements();
    
//    void simple_printf(const char *fmt, ...)
//{
//    va_list args;
//    va_start(args, fmt);
//
//    while (*fmt != '\0') {
//        if (*fmt == 'd') {
//            int i = va_arg(args, int);
//            std::cout << i << '\n';
//        } else if (*fmt == 's') {
//            char * s = va_arg(args, char*);
//            std::cout << s << '\n';
//        }
//        ++fmt;
//    }
//
//    va_end(args);
//}

    
};