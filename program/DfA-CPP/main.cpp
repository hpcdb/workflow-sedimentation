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
#include <vector>

/*
 * 
 */
int main(int argc, char** argv) {

    cout << "running main..." << endl;
    
    Dataflow dataflow = Dataflow("clothing");
    
    //set ideduplication
    Set ideduplication = Set("ideduplication");
    ideduplication.add_attribute("customerid", NUMERIC);
    ideduplication.add_attribute("country", TEXT);
    ideduplication.add_attribute("continent", TEXT);
    ideduplication.add_attribute("age", NUMERIC);
    ideduplication.add_attribute("gender", TEXT);
    ideduplication.add_attribute("children", NUMERIC);
    ideduplication.add_attribute("status", TEXT);
    dataflow.add_set(ideduplication);
    
    //set odeduplication
    Set odeduplication = Set("odeduplication");
    odeduplication.add_attribute("customerid", NUMERIC);
    odeduplication.add_attribute("country", TEXT);
    odeduplication.add_attribute("continent", TEXT);
    odeduplication.add_attribute("age", NUMERIC);
    odeduplication.add_attribute("gender", TEXT);
    odeduplication.add_attribute("children", NUMERIC);
    odeduplication.add_attribute("status", TEXT);
    dataflow.add_set(odeduplication);
    
    Transformation deduplication = Transformation("deduplication");
    dataflow.add_transformation(deduplication, ideduplication, odeduplication);
    
    dataflow.save();

    return 0;
}

