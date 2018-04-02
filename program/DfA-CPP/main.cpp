/* 
 * File:   main.cpp
 * Author: vitor
 *
 * Created on March 31, 2018, 10:55 AM
 */

#include "dataflow.h"
#include "task.h"

#include <cstdlib>
#include <string>

/*
 * 
 */
int main(int argc, char** argv) {

    cout << "running main..." << endl;
    
    Dataflow dataflow = Dataflow("clothing");
    dataflow.save();

    return 0;
}

